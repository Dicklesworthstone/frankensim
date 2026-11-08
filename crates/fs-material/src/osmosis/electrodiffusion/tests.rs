use super::*;
use crate::osmosis::{
    IdealDiluteSolution,
    donnan::{DonnanControls, DonnanProblem, VirialDonnan},
    nonideal::VirialSolutionModel,
};
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = 8.314_462_618_153_24 * T;

fn ideal(concentration: Vec<f64>) -> AqueousSolution {
    IdealDiluteSolution::new(T, 18e-6, 1.0, concentration)
        .unwrap()
        .into()
}
fn model(reference: f64, matrix: Vec<f64>, n: usize) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, 18e-6, reference, vec![1000.0; n], matrix).unwrap())
}
fn law(hydraulic: f64, valences: Vec<i32>, permeances: &[f64]) -> ElectrodiffusiveMembrane {
    ElectrodiffusiveMembrane::new(
        HydraulicPermeance::new(hydraulic).unwrap(),
        valences,
        permeances
            .iter()
            .map(|p| SolutePermeance::new(*p).unwrap())
            .collect(),
    )
    .unwrap()
}
fn problem<'a>(
    left: &'a AqueousSolution,
    right: &'a AqueousSolution,
) -> ElectrodiffusionProblem<'a> {
    ElectrodiffusionProblem {
        left,
        right,
        pressure_drop_pa: 0.0,
        trapped_particles_mol: [0.0; 2],
    }
}
fn psi(report: ElectrodiffusionReport) -> f64 {
    match report.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            ..
        } => dimensionless_left_minus_right,
        value => panic!("expected finite voltage, got {value:?}"),
    }
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "actual={actual:e}, expected={expected:e}, tolerance={absolute:e}"
    );
}

#[test]
fn g1_goldman_binary_voltage_and_flux_match_closed_form_with_fixed_charge_endpoints() {
    let (left, right) = (ideal(vec![1.2, 0.2]), ideal(vec![0.8, 0.8]));
    let membrane = law(1e-8, vec![1, -1], &[0.7, 0.2]);
    let mut flux = [0.0; 2];
    let report = membrane
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    // The general binary GHK voltage includes each endpoint concentration;
    // it need not assume that the MOBILE pool alone has zero bulk charge.
    let expected_psi = ((0.7_f64 * 0.8 + 0.2 * 0.2) / (0.7 * 1.2 + 0.2 * 0.8)).ln();
    near(psi(report), expected_psi, 2e-13);
    near(flux[0], -0.071_515_587_327_238_7, 2e-14);
    near(flux[1], flux[0], 1e-14);
    near(report.osmotic_pressure_drop_pa, -0.2 * RT, 2e-11);
    near(report.volume_flux_m_s, 1e-8 * 0.2 * RT, 2e-19);
    assert!(report.ionic_dissipation_w_m2.unwrap() > 0.0);
    assert!(report.hydraulic_dissipation_w_m2 > 0.0);
    assert!(report.potential_evaluations > 1);
}

#[test]
fn g1_neutral_salt_enters_exact_ion_free_bath_without_seeded_concentrations() {
    let (left, right) = (ideal(vec![100.0, 100.0]), ideal(vec![0.0, 0.0]));
    let membrane = law(1e-12, vec![1, -1], &[3e-6, 1e-6]);
    let mut flux = [0.0; 2];
    let report = membrane
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    near(psi(report), -(3.0_f64).ln(), 2e-12);
    near(flux[0], 1.5e-4 * 3.0_f64.ln(), 2e-16);
    near(flux[1], flux[0], 2e-16);
    assert!(report.volume_flux_m_s < 0.0);
    assert_eq!(report.ionic_dissipation_w_m2, None);
    assert_eq!(report.total_dissipation_w_m2, None);
    assert_eq!(right.solute_moles(), [0.0, 0.0]);
}

#[test]
fn g1_virial_donnan_partition_is_a_zero_flux_kinetic_equilibrium() {
    let shared = model(1.0, vec![0.8, 0.0, 0.0, 0.8], 2);
    let mut equilibrium = VirialDonnan::new(vec![1, -1], Arc::clone(&shared), 2e-13).unwrap();
    let (mut left_moles, mut right_moles) = ([0.0; 2], [0.0; 2]);
    let reference = equilibrium
        .equilibrate_into(
            DonnanProblem {
                temperature_k: T,
                volumes_m3: [1.0, 1.0],
                mobile_totals_mol: &[2.0, 1.0],
                fixed_charge_mol: [-1.0, 0.0],
                trapped_particles_mol: [0.1, 0.0],
            },
            DonnanControls::default(),
            &mut left_moles,
            &mut right_moles,
            || false,
        )
        .unwrap();
    let left = AqueousSolution::nonideal(Arc::clone(&shared), 1.0, left_moles.to_vec()).unwrap();
    let right = AqueousSolution::nonideal(shared, 1.0, right_moles.to_vec()).unwrap();
    let membrane = law(1e-8, vec![1, -1], &[0.7, 0.2]);
    let mut input = problem(&left, &right);
    input.trapped_particles_mol = [0.1, 0.0];
    input.pressure_drop_pa = reference.osmotic_pressure_pa[0] - reference.osmotic_pressure_pa[1];
    let mut flux = [0.0; 2];
    let report = membrane
        .flux_into(
            input,
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    near(flux[0], 0.0, 8e-14);
    near(flux[1], 0.0, 8e-14);
    near(report.volume_flux_m_s, 0.0, 3e-19);
    match reference.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            ..
        } => near(psi(report), dimensionless_left_minus_right, 4e-13),
        _ => panic!("ordinary salt has a finite equilibrium voltage"),
    }
    assert!(report.ionic_dissipation_w_m2.unwrap() >= 0.0);
}

#[test]
fn g3_multivalent_original_activity_flux_and_reversal_match_independent_oracle() {
    let shared = model(1.0, vec![0.8, 0.0, 0.0, 0.8], 2);
    let left = AqueousSolution::nonideal(
        Arc::clone(&shared),
        1.0,
        vec![0.736_039_179_709_984_6, 0.472_078_359_419_969_1],
    )
    .unwrap();
    let right = AqueousSolution::nonideal(
        shared,
        1.0,
        vec![0.763_960_820_290_015_4, 1.527_921_640_580_030_9],
    )
    .unwrap();
    let membrane = law(0.0, vec![2, -1], &[0.7, 0.2]);
    let (mut forward, mut reverse) = ([0.0; 2], [0.0; 2]);
    let a = membrane
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut forward,
            || false,
        )
        .unwrap();
    let b = membrane
        .flux_into(
            problem(&right, &left),
            ElectrodiffusionControls::default(),
            &mut reverse,
            || false,
        )
        .unwrap();
    near(psi(a), -0.141_687_308_484_333_3, 2e-13);
    near(forward[0], -0.180_195_898_549_922_57, 3e-14);
    near(forward[1], -0.360_391_797_099_845_26, 3e-14);
    near(
        a.ionic_dissipation_w_m2.unwrap(),
        1830.548_793_873_918_6,
        2e-9,
    );
    near(psi(a), -psi(b), 2e-13);
    for i in 0..2 {
        near(forward[i], -reverse[i], 2e-14);
    }
    near(
        a.osmotic_pressure_drop_pa,
        -b.osmotic_pressure_drop_pa,
        2e-10,
    );
}

#[test]
fn g0_counterion_limits_absent_channels_and_neutral_cosolutes_are_exact() {
    let (left, right) = (ideal(vec![1.0, 2.0]), ideal(vec![0.0, 1.0]));
    let membrane = law(1e-9, vec![1, 0], &[0.4, 0.2]);
    let mut flux = [9.0; 2];
    let report = membrane
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    assert_eq!(report.potential, DonnanPotential::NegativeInfinite);
    assert_eq!(flux[0], 0.0);
    near(flux[1], 0.2, 1e-15);
    near(
        report.ionic_dissipation_w_m2.unwrap(),
        RT * 0.2 * 2.0_f64.ln(),
        2e-12,
    );
    let reversed = membrane
        .flux_into(
            problem(&right, &left),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    assert_eq!(reversed.potential, DonnanPotential::PositiveInfinite);
    assert_eq!(flux[0], 0.0);
    near(flux[1], -0.2, 1e-15);
    let blocked = law(1e-9, vec![1, 0], &[0.0, 0.0]);
    let report = blocked
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    assert_eq!(report.potential, DonnanPotential::Indeterminate);
    assert_eq!(flux, [0.0; 2]);
    assert_eq!(report.ionic_dissipation_w_m2, Some(0.0));
    assert!(report.volume_flux_m_s < 0.0);
}

#[test]
fn g3_trace_currents_and_chemical_reference_changes_preserve_the_original_voltage() {
    let membrane = law(0.0, vec![1, -1], &[0.7, 0.2]);
    let (left, right) = (
        ideal(vec![1.2e-240, 0.2e-240]),
        ideal(vec![0.8e-240, 0.8e-240]),
    );
    let mut flux = [0.0; 2];
    let report = membrane
        .flux_into(
            problem(&left, &right),
            ElectrodiffusionControls::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    near(psi(report), 0.6_f64.ln(), 3e-12);
    near(flux[0] / 1e-240, -0.071_515_587_327_238_7, 4e-14);
    let evaluate = |reference: f64| {
        let shared = model(reference, vec![0.8, 0.0, 0.0, 0.8], 2);
        let left = AqueousSolution::nonideal(Arc::clone(&shared), 1.0, vec![1.2, 0.2]).unwrap();
        let right = AqueousSolution::nonideal(shared, 1.0, vec![0.8, 0.8]).unwrap();
        let mut flux = [0.0; 2];
        let report = membrane
            .flux_into(
                problem(&left, &right),
                ElectrodiffusionControls::default(),
                &mut flux,
                || false,
            )
            .unwrap();
        (flux, report)
    };
    assert_eq!(evaluate(1.0), evaluate(1e200));
}

#[test]
fn g4_failed_or_cancelled_trials_count_work_and_publish_no_partial_flux() {
    let (left, right) = (ideal(vec![1.2, 0.2]), ideal(vec![0.8, 0.8]));
    let membrane = law(1e-8, vec![1, -1], &[0.7, 0.2]);
    let sentinel = [17.0, -23.0];
    let mut output = sentinel;
    let mut work = 7;
    let controls = ElectrodiffusionControls {
        maximum_potential_evaluations: 1,
        ..Default::default()
    };
    assert_eq!(
        membrane.flux_into_counted(
            problem(&left, &right),
            controls,
            &mut output,
            &mut work,
            || false
        ),
        Err(ElectrodiffusionError::NonConvergence)
    );
    assert_eq!(work, 8);
    assert_eq!(output, sentinel);
    let mut polls = 0;
    let expected = membrane
        .flux_into_counted(
            problem(&left, &right),
            Default::default(),
            &mut output,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let completed_polls = polls;
    let mut polls = 0;
    output = sentinel;
    let before = work;
    assert_eq!(
        membrane.flux_into_counted(
            problem(&left, &right),
            Default::default(),
            &mut output,
            &mut work,
            || {
                polls += 1;
                polls == completed_polls
            }
        ),
        Err(ElectrodiffusionError::Cancelled)
    );
    assert_eq!(output, sentinel);
    assert_eq!(work - before, expected.potential_evaluations);
    let report = membrane
        .flux_into(
            problem(&left, &right),
            Default::default(),
            &mut output,
            || false,
        )
        .unwrap();
    assert_eq!(report, expected);
    let other =
        AqueousSolution::nonideal(model(1.0, vec![0.0; 4], 2), 1.0, vec![0.8, 0.8]).unwrap();
    output = sentinel;
    assert!(
        membrane
            .flux_into(
                problem(&left, &other),
                Default::default(),
                &mut output,
                || false
            )
            .is_err()
    );
    assert_eq!(output, sentinel);
}

#[test]
fn g1_channel_affinity_and_dissipation_use_identical_rounding_at_detailed_balance() {
    let p = 0.7;
    let left = 1.5_f64;
    let right = 0.5_f64;
    // A large cancelling excess/electric pair leaves a small log term. Both
    // diagnostics must use exactly the same compensated conjugate force.
    for electric in [-1e12 - left.ln() + right.ln(), -left.ln() + right.ln()] {
        let excess = if electric < -1e6 { 1e12 } else { 0.0 };
        let force = affinity(left, right, excess, electric).unwrap();
        let (flux, _) = channel_flux(p, excess + electric, left, right, Some(force)).unwrap();
        assert!(flux * force >= 0.0);
    }
}

#[test]
fn g1_analytic_response_matches_closed_binary_goldman_derivative_and_one_solve() {
    let (left, right) = (ideal(vec![3.0, 3.0]), ideal(vec![1.0, 1.0]));
    let membrane = law(1e-6, vec![1, -1], &[0.7, 0.2]);
    let mut flux = [0.0; 2];
    let mut work = 0;
    let response = membrane
        .flux_response_into_counted(
            problem(&left, &right),
            Default::default(),
            &mut flux,
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(work, response.report().potential_evaluations);
    let mut tangent = [0.0; 2];
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[1.0, 1.0],
                right_amounts_mol: &[-1.0, -1.0],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 13.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    // Differentiate the independent closed-form GHK voltage. The salt sum is
    // constant, so dJ=-P+P-/(P+-P-)*(cL+cR)*dpsi.
    let expected_psi = -0.5 / 1.3 - 0.5 / 2.3;
    let expected_flux = -(0.7 * 0.2 / 0.5) * 4.0 * expected_psi;
    near(
        report.dimensionless_potential_derivative.unwrap(),
        expected_psi,
        3e-13,
    );
    for value in tangent {
        near(value, expected_flux, 3e-13);
    }
    near(report.osmotic_pressure_drop_derivative_pa, 4.0 * RT, 2e-11);
    near(
        report.volume_flux_derivative_m_s,
        1e-6 * (13.0 - 4.0 * RT),
        2e-17,
    );
    near(report.charge_flux_derivative_mol_m2_s, 0.0, 2e-14);
    assert_eq!(work, response.report().potential_evaluations);
}

#[test]
fn g3_original_cross_virial_response_matches_full_inventory_volume_and_load_perturbation() {
    let mut matrix = vec![0.03; 36];
    for i in 0..6 {
        matrix[6 * i + i] = 0.4 + 0.08 * i as f64;
    }
    let shared = model(1.0, matrix, 6);
    let initial_left = [0.84, 0.14, 0.245, 0.28, 0.35, 0.0];
    let initial_right = [1.04, 1.04, 0.78, 0.91, 0.13, 0.0];
    let dl = [-0.08, 0.01, 0.02, 0.03, -0.025, 0.0];
    let dr = dl.map(|n| -n);
    let volumes = [0.7, 1.3];
    let dv = [-0.02, 0.02];
    let left =
        AqueousSolution::nonideal(Arc::clone(&shared), volumes[0], initial_left.to_vec()).unwrap();
    let right =
        AqueousSolution::nonideal(Arc::clone(&shared), volumes[1], initial_right.to_vec()).unwrap();
    let membrane = law(
        1e-6,
        vec![1, -1, 2, 0, -2, 1],
        &[0.7, 0.2, 0.4, 0.3, 0.0, 0.1],
    );
    let mut input = problem(&left, &right);
    input.trapped_particles_mol = [0.1, 0.03];
    let mut flux = [0.0; 6];
    let response = membrane
        .flux_response_into(input, Default::default(), &mut flux, || false)
        .unwrap();
    let mut tangent = [0.0; 6];
    let actual = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &dl,
                right_amounts_mol: &dr,
                volumes_m3: dv,
                pressure_drop_pa: 17.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    let evaluate = |epsilon: f64| {
        let left = AqueousSolution::nonideal(
            Arc::clone(&shared),
            volumes[0] + epsilon * dv[0],
            initial_left
                .iter()
                .zip(dl)
                .map(|(n, d)| n + epsilon * d)
                .collect(),
        )
        .unwrap();
        let right = AqueousSolution::nonideal(
            Arc::clone(&shared),
            volumes[1] + epsilon * dv[1],
            initial_right
                .iter()
                .zip(dr)
                .map(|(n, d)| n + epsilon * d)
                .collect(),
        )
        .unwrap();
        let mut input = problem(&left, &right);
        input.trapped_particles_mol = [0.1, 0.03];
        input.pressure_drop_pa = 17.0 * epsilon;
        let mut flux = [0.0; 6];
        let report = membrane
            .flux_into(input, Default::default(), &mut flux, || false)
            .unwrap();
        (flux, report)
    };
    let epsilon = 1e-5;
    let (plus, rplus) = evaluate(epsilon);
    let (minus, rminus) = evaluate(-epsilon);
    for i in 0..6 {
        near(tangent[i], (plus[i] - minus[i]) / (2.0 * epsilon), 3e-8);
    }
    assert_eq!(tangent[4], 0.0);
    assert_eq!(tangent[5], 0.0);
    near(
        actual.dimensionless_potential_derivative.unwrap(),
        (psi(rplus) - psi(rminus)) / (2.0 * epsilon),
        3e-8,
    );
    near(
        actual.volume_flux_derivative_m_s,
        (rplus.volume_flux_m_s - rminus.volume_flux_m_s) / (2.0 * epsilon),
        3e-11,
    );
    near(
        actual.osmotic_pressure_drop_derivative_pa,
        (rplus.osmotic_pressure_drop_pa - rminus.osmotic_pressure_drop_pa) / (2.0 * epsilon),
        3e-5,
    );
    near(actual.charge_flux_derivative_mol_m2_s, 0.0, 3e-14);
}

#[test]
fn g1_finite_voltage_derivative_includes_an_exact_absent_ion_without_log_floors() {
    let (left, right) = (ideal(vec![3.0, 3.0, 0.0]), ideal(vec![1.0, 1.0, 0.0]));
    let membrane = law(0.0, vec![1, -1, 2], &[0.7, 0.2, 0.3]);
    let mut flux = [0.0; 3];
    let response = membrane
        .flux_response_into(
            problem(&left, &right),
            Default::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    assert_eq!(flux[2], 0.0);
    let mut tangent = [0.0; 3];
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[0.0, 0.0, 1.0],
                right_amounts_mol: &[0.0; 3],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    let a = 2.0 * psi(*response.report());
    let donor = 0.3 * (-a) / (-a).exp_m1();
    near(tangent[2], donor, 2e-14);
    near(tangent[0] - tangent[1] + 2.0 * tangent[2], 0.0, 2e-13);
    assert!(report.dimensionless_potential_derivative.unwrap() < 0.0);
    // One-sided salt endpoints likewise have finite original stream derivatives.
    let (left, right) = (ideal(vec![1.0, 1.0]), ideal(vec![0.0, 0.0]));
    let membrane = law(0.0, vec![1, -1], &[0.7, 0.2]);
    let response = membrane
        .flux_response_into(
            problem(&left, &right),
            Default::default(),
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[-1.0, -1.0],
                right_amounts_mol: &[1.0, 1.0],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    near(
        report.dimensionless_potential_derivative.unwrap(),
        0.5 / 0.2 + 0.5 / 0.7,
        3e-12,
    );
}

#[test]
fn g0_singular_voltage_allows_fixed_support_water_neutral_tangents_only() {
    let (left, right) = (ideal(vec![1.0, 2.0]), ideal(vec![0.0, 1.0]));
    let membrane = law(1e-6, vec![1, 0], &[0.7, 0.2]);
    let response = membrane
        .flux_response_into(
            problem(&left, &right),
            Default::default(),
            &mut [0.0; 2],
            || false,
        )
        .unwrap();
    let mut tangent = [0.0; 2];
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[0.0, 1.0],
                right_amounts_mol: &[0.0, -1.0],
                volumes_m3: [-0.1, 0.1],
                pressure_drop_pa: 7.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    assert_eq!(report.dimensionless_potential_derivative, None);
    assert_eq!(tangent[0], 0.0);
    near(tangent[1], 0.46, 3e-15);
    near(
        report.volume_flux_derivative_m_s,
        1e-6 * (7.0 - 2.4 * RT),
        2e-17,
    );
    tangent = [17.0, -23.0];
    assert_eq!(
        response.directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[-1.0, 0.0],
                right_amounts_mol: &[1.0, 0.0],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut tangent,
            || false
        ),
        Err(ElectrodiffusionError::SingularDerivative)
    );
    assert_eq!(tangent, [17.0, -23.0]);
    let (left, right) = (ideal(vec![0.0, 0.0, 2.0]), ideal(vec![0.0, 0.0, 1.0]));
    let membrane = law(0.0, vec![1, -1, 0], &[0.7, 0.2, 0.4]);
    let response = membrane
        .flux_response_into(
            problem(&left, &right),
            Default::default(),
            &mut [0.0; 3],
            || false,
        )
        .unwrap();
    let mut tangent = [0.0; 3];
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[0.0, 0.0, 0.1],
                right_amounts_mol: &[0.0, 0.0, -0.1],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    assert_eq!(report.dimensionless_potential_derivative, None);
    near(tangent[2], 0.08, 2e-15);
}

#[test]
fn g3_charge_scaling_ignores_large_neutral_rates_and_harmless_absent_donor_overflow() {
    let (left, right) = (
        ideal(vec![3e-200, 3e-200, 1e100, 0.0]),
        ideal(vec![1e-200, 1e-200, 1e100, 0.0]),
    );
    let membrane = law(
        0.0,
        vec![1, -1, 0, 2_000_000_000],
        &[0.7, 0.2, 1e100, 1e308],
    );
    let response = membrane
        .flux_response_into(
            problem(&left, &right),
            Default::default(),
            &mut [0.0; 4],
            || false,
        )
        .unwrap();
    let mut tangent = [0.0; 4];
    let report = response
        .directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[1e-200, 1e-200, 1e100, 0.0],
                right_amounts_mol: &[-1e-200, -1e-200, -1e100, 0.0],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut tangent,
            || false,
        )
        .unwrap();
    near(
        report.dimensionless_potential_derivative.unwrap(),
        -0.5 / 1.3 - 0.5 / 2.3,
        3e-12,
    );
    near(tangent[2] / 1e200, 2.0, 2e-14);
    assert_eq!(tangent[3], 0.0);
    assert!(matches!(
        response.directional_derivative_into(
            ElectrodiffusionDirection {
                left_amounts_mol: &[0.0, 0.0, 0.0, 1.0],
                right_amounts_mol: &[0.0; 4],
                volumes_m3: [0.0; 2],
                pressure_drop_pa: 0.0,
            },
            &mut tangent,
            || false
        ),
        Err(ElectrodiffusionError::SingularDerivative)
    ));
}

#[test]
fn g4_response_and_tangent_cancel_without_publishing_and_reuse_original_work() {
    let (left, right) = (ideal(vec![3.0, 3.0]), ideal(vec![1.0, 1.0]));
    let membrane = law(1e-6, vec![1, -1], &[0.7, 0.2]);
    let mut output = [17.0, -23.0];
    let mut work = 0;
    let mut polls = 0;
    let response = membrane
        .flux_response_into_counted(
            problem(&left, &right),
            Default::default(),
            &mut output,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let total_polls = polls;
    let expected_work = work;
    let direction = ElectrodiffusionDirection {
        left_amounts_mol: &[1.0, 1.0],
        right_amounts_mol: &[-1.0, -1.0],
        volumes_m3: [0.0; 2],
        pressure_drop_pa: 0.0,
    };
    let mut tangent = [17.0, -23.0];
    let mut derivative_polls = 0;
    let expected = response
        .directional_derivative_into(direction, &mut tangent, || {
            derivative_polls += 1;
            false
        })
        .unwrap();
    tangent = [17.0, -23.0];
    let mut polls = 0;
    assert_eq!(
        response.directional_derivative_into(direction, &mut tangent, || {
            polls += 1;
            polls == derivative_polls
        }),
        Err(ElectrodiffusionError::Cancelled)
    );
    assert_eq!(tangent, [17.0, -23.0]);
    assert_eq!(
        response
            .directional_derivative_into(direction, &mut tangent, || false)
            .unwrap(),
        expected
    );
    output = [17.0, -23.0];
    let mut polls = 0;
    let mut work = 0;
    assert!(matches!(
        membrane.flux_response_into_counted(
            problem(&left, &right),
            Default::default(),
            &mut output,
            &mut work,
            || {
                polls += 1;
                polls == total_polls
            }
        ),
        Err(ElectrodiffusionError::Cancelled)
    ));
    assert_eq!(output, [17.0, -23.0]);
    assert_eq!(work, expected_work);
}
