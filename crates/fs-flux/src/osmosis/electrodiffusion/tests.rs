//! Independent finite-rate electrochemical, mechanical and temporal regressions.

use super::*;
use crate::osmosis::{OsmoticCompartment, OsmoticPressureVolume};
use fs_material::{
    hyper::{
        Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric,
        discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings},
    },
    osmosis::{
        AqueousSolution, HydraulicPermeance, IdealDiluteSolution, SolutePermeance,
        donnan::DonnanPotential,
        electrodiffusion::{
            ElectrodiffusionControls, ElectrodiffusionReport, ElectrodiffusiveMembrane,
        },
        nonideal::VirialSolutionModel,
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = R * T;

#[test]
fn small_ideal_salt_transfer_uses_its_linear_root_without_bath_scale_bisection() {
    let mut pair = Case {
        amounts: [[2.0, 2.0], [1.0, 1.0]],
        interaction: None,
        permeances: [1.0, 1.0],
        hydraulic: 0.0,
        fixed: [0.0; 2],
        trapped: [0.0; 2],
        ..Case::default()
    }
    .build();
    let dt = 1e-6;
    let report = pair.step(dt, controls(), || false).unwrap();
    // Equal ion permeances and neutral baths give psi=0 and the independent
    // finite-reservoir BE equation q=dt*(1-2q). Its root is near q=0, while
    // the physically allowed bracket contains nearly two whole moles.
    let expected = dt / (1.0 + 2.0 * dt);
    for transfer in report.solute_transfer_mol {
        near(transfer, expected, 3e-15);
    }
    assert_eq!(report.solvent_transfer_m3, 0.0);
    assert_eq!(psi(&report.endpoint), 0.0);
    assert!(report.salt_evaluations <= 4);
}

#[derive(Clone, Copy)]
struct Case {
    amounts: [[f64; 2]; 2],
    valences: [i32; 2],
    interaction: Option<f64>,
    permeances: [f64; 2],
    hydraulic: f64,
    fixed: [f64; 2],
    trapped: [f64; 2],
    gauge: f64,
    initial_electrical: Option<ElectrodiffusionControls>,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            amounts: [[1.2, 0.2], [0.8, 0.8]],
            valences: [1, -1],
            interaction: Some(0.8),
            permeances: [0.7, 0.2],
            hydraulic: 1.0 / RT,
            fixed: [-1.0, 0.0],
            trapped: [0.1, 0.0],
            gauge: 0.0,
            initial_electrical: None,
        }
    }
}

fn pressure(gauge: f64) -> OsmoticPressureVolume {
    let material = Hyperelastic::new(
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 1.0,
                alpha: 2.0,
            }],
            volumetric: Volumetric::SimoTaylor { kappa: RT },
        },
        1.0,
    )
    .unwrap();
    let point = DiscreteBulkPoint::try_new(
        &material,
        1.0,
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.01,
            maximum_volume_ratio: 3.0,
            maximum_iterations: 96,
            absolute_pressure_tolerance_pa: 1e-10,
            relative_pressure_tolerance: 1e-13,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-13,
        },
    )
    .unwrap();
    OsmoticPressureVolume::from_bulk_material(1.0, 0.0, gauge, point).unwrap()
}

impl Case {
    fn build(self) -> IonicOsmoticPair {
        let model = self.interaction.map(|b| {
            Arc::new(
                VirialSolutionModel::new(T, 18e-6, 1.0, vec![128.0; 2], vec![b, 0.0, 0.0, b])
                    .unwrap(),
            )
        });
        let cells = self.amounts.map(|amounts| {
            let solution = match &model {
                Some(model) => {
                    AqueousSolution::nonideal(Arc::clone(model), 1.0, amounts.to_vec()).unwrap()
                }
                None => IdealDiluteSolution::new(T, 18e-6, 1.0, amounts.to_vec())
                    .unwrap()
                    .into(),
            };
            OsmoticCompartment::new(solution, pressure(self.gauge), 0.05).unwrap()
        });
        let membrane = ElectrodiffusiveMembrane::new(
            HydraulicPermeance::new(self.hydraulic).unwrap(),
            self.valences.to_vec(),
            self.permeances
                .into_iter()
                .map(|p| SolutePermeance::new(p).unwrap())
                .collect(),
        )
        .unwrap();
        match self.initial_electrical {
            Some(electrical) => IonicOsmoticPair::new_with_electrical_controls(
                cells,
                membrane,
                self.fixed,
                self.trapped,
                1.0,
                1e-11,
                electrical,
                1 << 20,
            ),
            None => IonicOsmoticPair::new(
                cells,
                membrane,
                self.fixed,
                self.trapped,
                1.0,
                1e-11,
                1 << 20,
            ),
        }
        .unwrap()
    }
}

fn controls() -> IonicStepControls {
    IonicStepControls {
        absolute_water_tolerance_m3: 1e-11,
        absolute_amount_tolerance_mol: 1e-11,
        relative_transfer_tolerance: 1e-11,
        free_energy_tolerance_j: 1e-7,
        ..IonicStepControls::default()
    }
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

fn psi(report: &ElectrodiffusionReport) -> f64 {
    match report.potential {
        DonnanPotential::Finite {
            dimensionless_left_minus_right,
            ..
        } => dimensionless_left_minus_right,
        other => panic!("expected a finite voltage, got {other:?}"),
    }
}

fn snapshot(pair: &IonicOsmoticPair) -> (Vec<u64>, ElectrodiffusionReport) {
    let mut bits = vec![pair.time_s().to_bits()];
    for cell in pair.compartments() {
        bits.push(cell.solution().solvent_volume_m3().to_bits());
        bits.extend(cell.solution().solute_moles().iter().map(|n| n.to_bits()));
    }
    bits.extend(pair.mobile_totals_mol().iter().map(|n| n.to_bits()));
    bits.extend(pair.fixed_charge_mol().iter().map(|n| n.to_bits()));
    bits.extend(pair.trapped_particles_mol().iter().map(|n| n.to_bits()));
    (bits, *pair.endpoint())
}

// Original analytic potentials, without calls to the production energy owner.
// The reference is c_ref=1; conserved reference terms cancel from every change.
fn potential(pair: &IonicOsmoticPair, b: f64, trapped: [f64; 2]) -> f64 {
    let mut value = 0.0;
    for (side, cell) in pair.compartments().iter().enumerate() {
        let v = cell.solution().solvent_volume_m3();
        for &n in cell.solution().solute_moles() {
            if n > 0.0 {
                value += n * ((n / v).ln() - 1.0) + 0.5 * b * n * n / v;
            }
        }
        if trapped[side] > 0.0 {
            value += trapped[side] * ((trapped[side] / v).ln() - 1.0);
        }
        value += 0.25 * (v * v - 1.0 - 2.0 * v.ln());
    }
    RT * value
}

#[test]
fn nonideal_backward_euler_matches_independent_coupled_oracle_and_pressure_gauge() {
    // Independent scalar bisections solved current, salt BE and water BE from
    // GHK, the virial potential, and original Simo--Taylor U and U'. Salt entry
    // reverses the initial water direction; an instantaneous Donnan partition
    // would instead end at x=-0.08159122061597, t=0.35913361541131.
    let case = Case::default();
    let mut pair = case.build();
    assert_eq!(
        pair.compartments()[0].solution().solute_moles(),
        &[1.2, 0.2]
    );
    assert!(pair.endpoint().volume_flux_m_s > 0.0);
    let before = potential(&pair, 0.8, case.trapped);
    let report = pair.step(0.2, controls(), || false).unwrap();
    let x = report.solvent_transfer_m3;
    near(x, -0.006_246_802_818_186_015, 2e-9);
    let expected = [
        [1.215_132_808_256_861_5, 0.215_132_808_256_861_6],
        [0.784_867_191_743_138_5, 0.784_867_191_743_138_5],
    ];
    for (side, cell) in pair.compartments().iter().enumerate() {
        for (i, &n) in cell.solution().solute_moles().iter().enumerate() {
            near(n, expected[side][i], 2e-9);
        }
    }
    near(psi(&report.endpoint), -0.870_785_996_090_213_3, 2e-9);
    near(report.endpoint.volume_flux_m_s, x / 0.2, 2e-9);
    near(
        report.endpoint.driving_pressure_pa,
        -77.427_778_793_409_5,
        2e-5,
    );
    for transfer in report.solute_transfer_mol {
        near(transfer / 0.2, -0.075_664_041_284_308, 2e-8);
    }
    near(report.water_balance_residual_m3, 0.0, 1e-13);
    near(report.water_root_residual_m3, 0.0, 1e-9);
    for residual in report.solute_balance_residual_mol {
        near(residual, 0.0, 1e-13);
    }
    for residual in report.solute_root_residual_mol {
        near(residual, 0.0, 1e-9);
    }
    near(report.endpoint.charge_flux_residual_mol_m2_s, 0.0, 2e-13);
    let exact_work = -0.5 * RT * (x * x - (-(x * x)).ln_1p());
    near(report.pressure_work_into_solution_j, exact_work, 2e-9);
    near(
        report.pressure_work_into_solution_j,
        -0.096_736_157_092_629_2,
        2e-6,
    );
    let change = potential(&pair, 0.8, case.trapped) - before;
    near(report.total_potential_change_j, change, 2e-8);
    near(change, -40.390_753_045_838_54, 2e-6);
    let dissipation = report.endpoint.total_dissipation_w_m2.unwrap();
    near(dissipation, 191.622_751_800_718_36, 2e-5);
    assert!(change <= -0.2 * dissipation + 1e-7);

    let mut shifted = Case {
        gauge: 1e24,
        ..case
    }
    .build();
    let translated = shifted.step(0.2, controls(), || false).unwrap();
    assert_eq!(snapshot(&pair), snapshot(&shifted));
    assert_eq!(
        report.pressure_work_into_solution_j.to_bits(),
        translated.pressure_work_into_solution_j.to_bits()
    );
    assert_eq!(
        report.total_potential_change_j.to_bits(),
        translated.total_potential_change_j.to_bits()
    );
}

#[test]
fn ideal_ghk_retains_finite_kinetics_and_approaches_the_donnan_limit() {
    // At fixed unit volumes, exact ideal open voltage is
    // ln[(P+ cR,+ + P- cL,-)/(P+ cL,+ + P- cR,-)]. With fixed charge -1,
    // the independent salt equation is t-.2 + dt*.14*B(psi)*(3t-1)/(.9+.5t)=0.
    for (scale, expected_t) in [
        (1e-6, 0.200_000_014_303_115_98),
        (1.0, 0.212_934_353_940_940_18),
        (1e3, 0.332_122_867_025_791_4),
    ] {
        let mut pair = Case {
            interaction: None,
            hydraulic: 0.0,
            permeances: [0.7 * scale, 0.2 * scale],
            trapped: [0.0; 2],
            initial_electrical: (scale > 1.0).then_some(ElectrodiffusionControls {
                absolute_current_tolerance_mol_m2_s: 1e-10,
                ..ElectrodiffusionControls::default()
            }),
            ..Case::default()
        }
        .build();
        near(psi(pair.endpoint()), 0.6_f64.ln(), 2e-11);
        let physical = IonicStepControls {
            relative_transfer_tolerance: 1e-10,
            electrical: ElectrodiffusionControls {
                absolute_current_tolerance_mol_m2_s: 1e-10,
                ..ElectrodiffusionControls::default()
            },
            ..controls()
        };
        let report = pair.step(0.2, physical, || false).unwrap();
        let t = pair.compartments()[0].solution().solute_moles()[1];
        near(t, expected_t, 2e-9);
        assert_eq!(report.solvent_transfer_m3, 0.0);
        assert_eq!(report.pressure_work_into_solution_j, 0.0);
        let expected_psi = ((0.7 - 0.5 * t) / (0.9 + 0.5 * t)).ln();
        near(psi(&report.endpoint), expected_psi, 2e-10);
        let bernoulli = expected_psi / expected_psi.exp_m1();
        let flux = scale * 0.14 * bernoulli * (3.0 * t - 1.0) / (0.9 + 0.5 * t);
        for transfer in report.solute_transfer_mol {
            near(transfer, 0.2 * flux, 2e-9);
        }
        if scale == 1.0 {
            assert!((t - 1.0 / 3.0).abs() > 0.1);
        } else if scale > 1.0 {
            assert!((t - 1.0 / 3.0).abs() < 0.002);
        }
    }
}

#[test]
fn blocked_ion_channels_and_absent_coions_do_not_invent_salt() {
    for permeances in [[0.0, 0.0], [0.7, 0.0]] {
        let mut pair = Case {
            permeances,
            ..Case::default()
        }
        .build();
        let report = pair.step(0.2, controls(), || false).unwrap();
        near(
            report.solvent_transfer_m3,
            0.001_624_926_373_615_831_5,
            2e-9,
        );
        for (side, cell) in pair.compartments().iter().enumerate() {
            for (i, &n) in cell.solution().solute_moles().iter().enumerate() {
                near(n, Case::default().amounts[side][i], 1e-10);
            }
        }
        for transfer in report.solute_transfer_mol {
            near(transfer, 0.0, 1e-10);
        }
        if permeances[0] == 0.0 {
            assert_eq!(report.endpoint.potential, DonnanPotential::Indeterminate);
        } else {
            let left = pair.compartments()[0].solution();
            let right = pair.compartments()[1].solution();
            let cl = left.solute_moles()[0] / left.solvent_volume_m3();
            let cr = right.solute_moles()[0] / right.solvent_volume_m3();
            near(
                psi(&report.endpoint),
                -(cl / cr).ln() - 0.8 * (cl - cr),
                2e-9,
            );
        }
    }
    let mut pair = Case {
        amounts: [[1.0, 0.0], [0.0, 0.0]],
        trapped: [0.0; 2],
        ..Case::default()
    }
    .build();
    let report = pair.step(0.2, controls(), || false).unwrap();
    near(report.solvent_transfer_m3, -0.164_216_010_241_167_44, 2e-9);
    assert_eq!(
        pair.compartments()[0].solution().solute_moles(),
        &[1.0, 0.0]
    );
    assert_eq!(
        pair.compartments()[1].solution().solute_moles(),
        &[0.0, 0.0]
    );
    assert_eq!(report.solute_transfer_mol, [0.0; 2]);
    assert_eq!(report.endpoint.potential, DonnanPotential::NegativeInfinite);
    assert!(report.total_potential_change_j < 0.0);
}

#[test]
fn supported_coupled_endpoint_can_require_dilution_before_any_salt_advance() {
    // Initially the left cation is exactly at its declared cap. At fixed water,
    // the independent ion BE root would transfer -0.000837184926569821 mol,
    // which is outside that cap for EVERY positive step duration. The complete
    // coupled endpoint instead dilutes first and remains strictly supported.
    // Halving time without searching the coupled domain cannot solve this case.
    let caps = [1.2, 1.0];
    let model = Arc::new(
        VirialSolutionModel::new(T, 18e-6, 1.0, caps.to_vec(), vec![0.8, 0.0, 0.0, 0.8]).unwrap(),
    );
    let amounts = [[1.2, 0.2], [0.8, 0.8]];
    let cells = [0, 1].map(|side| {
        OsmoticCompartment::new(
            AqueousSolution::nonideal(Arc::clone(&model), 1.0, amounts[side].to_vec()).unwrap(),
            OsmoticPressureVolume::new(1.0, [0.0, RT][side], 0.0).unwrap(),
            0.01,
        )
        .unwrap()
    });
    let membrane = ElectrodiffusiveMembrane::new(
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        vec![1, -1],
        vec![
            SolutePermeance::new(0.7).unwrap(),
            SolutePermeance::new(0.2).unwrap(),
        ],
    )
    .unwrap();
    let mut pair =
        IonicOsmoticPair::new(cells, membrane, [-1.0, 0.0], [0.0; 2], 1.0, 1e-11, 1 << 20).unwrap();
    let report = pair.step(0.01, controls(), || false).unwrap();
    near(report.solvent_transfer_m3, -0.008_417_969_072_558_327, 2e-9);
    for transfer in report.solute_transfer_mol {
        near(transfer, -0.000_880_660_644_880_319, 2e-9);
    }
    near(psi(&report.endpoint), -0.828_210_014_185_601_7, 2e-9);
    for cell in pair.compartments() {
        for (i, &n) in cell.solution().solute_moles().iter().enumerate() {
            assert!(n / cell.solution().solvent_volume_m3() <= caps[i]);
        }
    }
    near(
        pair.compartments()[0].solution().solute_moles()[0],
        1.200_880_660_644_880_3,
        2e-9,
    );
    // The imposed right pressure supplies work while the shared potential falls.
    near(
        report.pressure_work_into_solution_j,
        20.867_783_607_393_96,
        2e-6,
    );
    near(
        report.total_potential_change_j,
        -20.503_347_687_013_104,
        2e-6,
    );
    near(report.water_root_residual_m3, 0.0, 1e-9);
    for residual in report.solute_root_residual_mol {
        near(residual, 0.0, 1e-9);
    }
}

#[test]
fn unequal_valence_uses_neutral_stoichiometry_and_both_original_ion_equations() {
    // Independent current and salt bisections for a +2/-1 electrolyte.
    // Neutral transfer is [q,2q]; equal species transfers would create charge.
    let case = Case {
        amounts: [[0.7, 0.4], [0.8, 1.6]],
        valences: [2, -1],
        hydraulic: 0.0,
        trapped: [0.0; 2],
        ..Case::default()
    };
    let mut pair = case.build();
    let before = potential(&pair, 0.8, case.trapped);
    let report = pair.step(0.2, controls(), || false).unwrap();
    near(
        report.solute_transfer_mol[0],
        -0.036_039_179_709_984_55,
        2e-9,
    );
    near(
        report.solute_transfer_mol[1],
        -0.072_078_359_419_969_1,
        2e-9,
    );
    near(
        report.solute_transfer_mol[1],
        2.0 * report.solute_transfer_mol[0],
        1e-13,
    );
    near(psi(&report.endpoint), -0.141_687_308_484_333_3, 2e-9);
    for (side, cell) in pair.compartments().iter().enumerate() {
        let n = cell.solution().solute_moles();
        near(2.0 * n[0] - n[1] + case.fixed[side], 0.0, 1e-13);
    }
    for residual in report.solute_root_residual_mol {
        near(residual, 0.0, 1e-9);
    }
    near(report.endpoint.charge_flux_residual_mol_m2_s, 0.0, 2e-13);
    assert_eq!(report.solvent_transfer_m3, 0.0);
    near(
        report.total_potential_change_j,
        -401.830_379_365_856_06,
        2e-5,
    );
    near(
        report.total_potential_change_j,
        potential(&pair, 0.8, case.trapped) - before,
        2e-8,
    );
    // Reversing the declared species order must also reverse stoichiometric
    // transfer, without changing the physical voltage or the accepted state.
    let mut reversed = Case {
        amounts: case.amounts.map(|n| [n[1], n[0]]),
        valences: [-1, 2],
        permeances: [0.2, 0.7],
        ..case
    }
    .build();
    let reordered = reversed.step(0.2, controls(), || false).unwrap();
    near(psi(&reordered.endpoint), psi(&report.endpoint), 2e-9);
    for i in 0..2 {
        near(
            reordered.solute_transfer_mol[i],
            report.solute_transfer_mol[1 - i],
            2e-9,
        );
    }
    near(
        reordered.total_potential_change_j,
        report.total_potential_change_j,
        2e-7,
    );
}

#[test]
fn represented_clock_interval_drives_water_and_both_ion_equations() {
    let initial = Case::default().build();
    let mut ordinary = initial.clone();
    let expected = ordinary.step(2.0, controls(), || false).unwrap();
    let mut translated = initial;
    let clock = 4_503_599_627_370_496.0; // 2^52 seconds: adding 1.5 rounds to 2.
    translated.time_s = clock;
    let report = translated.step(1.5, controls(), || false).unwrap();
    assert_eq!(report.interval_s, [clock, clock + 2.0]);
    assert_eq!(translated.time_s(), clock + 2.0);
    for side in 0..2 {
        assert_eq!(
            translated.compartments()[side].solution(),
            ordinary.compartments()[side].solution(),
        );
    }
    assert_eq!(report.endpoint, expected.endpoint);
    assert_eq!(report.solvent_transfer_m3, expected.solvent_transfer_m3);
    assert_eq!(report.solute_transfer_mol, expected.solute_transfer_mol);
    assert_eq!(
        report.pressure_work_into_solution_j,
        expected.pressure_work_into_solution_j
    );
    assert_eq!(
        report.total_potential_change_j,
        expected.total_potential_change_j
    );
}

#[test]
fn adaptive_late_cancellation_restores_the_full_state_and_retries_exactly() {
    let physical = IonicAdvanceControls {
        maximum_step_s: 0.1,
        relative_tolerance: 5e-3,
        absolute_amount_tolerance_mol: 1e-8,
        absolute_water_tolerance_m3: 1e-8,
        step: controls(),
        ..IonicAdvanceControls::default()
    };
    let mut successful = Case::default().build();
    let mut successful_workspace = successful.adaptive_workspace(1 << 20).unwrap();
    let mut successful_polls = 0usize;
    let accepted = successful
        .advance(0.3, physical, &mut successful_workspace, || {
            successful_polls += 1;
            false
        })
        .unwrap();
    assert!(accepted.accepted_steps >= 3);
    assert!(accepted.maximum_scaled_error <= 1.0);
    assert!(accepted.endpoint.interval_s[0] > 0.0);
    assert_eq!(successful.time_s(), 0.3);
    assert!(successful_polls > 10);

    let mut pair = Case::default().build();
    let original = snapshot(&pair);
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let mut polls = 0usize;
    let mut cancellation_reached = false;
    let refused = pair.advance(0.3, physical, &mut workspace, || {
        polls += 1;
        cancellation_reached = polls >= successful_polls - 1;
        cancellation_reached
    });
    assert!(cancellation_reached);
    assert!(matches!(
        refused,
        Err(IonicOsmoticError::Network(OsmoticNetworkError::Cancelled))
    ));
    assert_eq!(snapshot(&pair), original);

    let retry = pair
        .advance(0.3, physical, &mut workspace, || false)
        .unwrap();
    assert_eq!(snapshot(&pair), snapshot(&successful));
    assert_eq!(retry.accepted_steps, accepted.accepted_steps);
    assert_eq!(retry.rejected_steps, accepted.rejected_steps);
    assert_eq!(
        retry.pressure_work_into_solution_j.to_bits(),
        accepted.pressure_work_into_solution_j.to_bits()
    );
}

#[derive(Clone)]
struct MixtureCase {
    amounts: [Vec<f64>; 2],
    valences: Vec<i32>,
    permeances: Vec<f64>,
    matrix: Option<Vec<f64>>,
    hydraulic: f64,
    fixed: [f64; 2],
    trapped: [f64; 2],
}

impl MixtureCase {
    fn sodium_potassium_chloride() -> Self {
        Self {
            amounts: [vec![1.0, 0.2, 0.2], vec![0.2, 0.6, 0.8]],
            valences: vec![1, 1, -1],
            permeances: vec![0.3, 0.8, 0.2],
            matrix: Some(vec![0.8, 0.15, 0.1, 0.15, 0.5, 0.12, 0.1, 0.12, 0.6]),
            hydraulic: 1.0 / RT,
            fixed: [-1.0, 0.0],
            trapped: [0.1, 0.0],
        }
    }

    fn build(&self) -> IonicOsmoticPair {
        let count = self.valences.len();
        let model = self.matrix.as_ref().map(|matrix| {
            Arc::new(
                VirialSolutionModel::new(T, 18e-6, 1.0, vec![128.0; count], matrix.clone())
                    .unwrap(),
            )
        });
        let cells = [0, 1].map(|side| {
            let solution = match &model {
                Some(model) => {
                    AqueousSolution::nonideal(Arc::clone(model), 1.0, self.amounts[side].clone())
                        .unwrap()
                }
                None => IdealDiluteSolution::new(T, 18e-6, 1.0, self.amounts[side].clone())
                    .unwrap()
                    .into(),
            };
            OsmoticCompartment::new(solution, pressure(0.0), 0.05).unwrap()
        });
        let membrane = ElectrodiffusiveMembrane::new(
            HydraulicPermeance::new(self.hydraulic).unwrap(),
            self.valences.clone(),
            self.permeances
                .iter()
                .map(|&p| SolutePermeance::new(p).unwrap())
                .collect(),
        )
        .unwrap();
        IonicOsmoticPair::new(
            cells,
            membrane,
            self.fixed,
            self.trapped,
            1.0,
            1e-11,
            1 << 20,
        )
        .unwrap()
    }

    fn permuted(&self, order: &[usize]) -> Self {
        let n = self.valences.len();
        assert_eq!(order.len(), n);
        let matrix = self.matrix.as_ref().map(|original| {
            let mut permuted = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    permuted[i * n + j] = original[order[i] * n + order[j]];
                }
            }
            permuted
        });
        Self {
            amounts: self
                .amounts
                .each_ref()
                .map(|ns| order.iter().map(|&i| ns[i]).collect()),
            valences: order.iter().map(|&i| self.valences[i]).collect(),
            permeances: order.iter().map(|&i| self.permeances[i]).collect(),
            matrix,
            hydraulic: self.hydraulic,
            fixed: self.fixed,
            trapped: self.trapped,
        }
    }
}

fn mixture_potential(pair: &IonicOsmoticPair, case: &MixtureCase) -> f64 {
    let mut value = 0.0;
    for (side, cell) in pair.compartments().iter().enumerate() {
        let n = cell.solution().solute_moles();
        let v = cell.solution().solvent_volume_m3();
        for &amount in n {
            if amount > 0.0 {
                value += amount * ((amount / v).ln() - 1.0);
            }
        }
        if let Some(matrix) = &case.matrix {
            for (i, &a) in n.iter().enumerate() {
                for (j, &b) in n.iter().enumerate() {
                    value += 0.5 * matrix[i * n.len() + j] * a * b / v;
                }
            }
        }
        if case.trapped[side] > 0.0 {
            value += case.trapped[side] * ((case.trapped[side] / v).ln() - 1.0);
        }
        value += 0.25 * (v * v - 1.0 - 2.0 * v.ln());
    }
    RT * value
}

fn independent_bernoulli(a: f64) -> f64 {
    if a.abs() < 1e-4 {
        let a2 = a * a;
        1.0 - 0.5 * a + a2 / 12.0 - a2 * a2 / 720.0
    } else {
        a / a.exp_m1()
    }
}

fn independent_mixture_fluxes(
    pair: &IonicOsmoticPair,
    case: &MixtureCase,
    voltage: f64,
) -> Vec<f64> {
    let c: [Vec<f64>; 2] = pair.compartments().each_ref().map(|cell| {
        cell.solution()
            .solute_moles()
            .iter()
            .map(|n| n / cell.solution().solvent_volume_m3())
            .collect()
    });
    let n = case.valences.len();
    (0..n)
        .map(|i| {
            let excess = case.matrix.as_ref().map_or(0.0, |b| {
                (0..n).map(|j| b[i * n + j] * (c[0][j] - c[1][j])).sum()
            });
            let a = excess + f64::from(case.valences[i]) * voltage;
            case.permeances[i]
                * (independent_bernoulli(-a) * c[0][i] - independent_bernoulli(a) * c[1][i])
        })
        .collect()
}

fn assert_mixture_conserved(pair: &IonicOsmoticPair, case: &MixtureCase) {
    assert_eq!(pair.species_count(), case.valences.len());
    assert_eq!(pair.mobile_totals_mol().len(), case.valences.len());
    for i in 0..case.valences.len() {
        let total: f64 = pair
            .compartments()
            .iter()
            .map(|c| c.solution().solute_moles()[i])
            .sum();
        near(total, case.amounts[0][i] + case.amounts[1][i], 5e-13);
    }
    for (side, cell) in pair.compartments().iter().enumerate() {
        let mobile_charge: f64 = cell
            .solution()
            .solute_moles()
            .iter()
            .zip(&case.valences)
            .map(|(&n, &z)| n * f64::from(z))
            .sum();
        near(mobile_charge + case.fixed[side], 0.0, 2e-11);
    }
    near(
        pair.compartments()
            .iter()
            .map(|c| c.solution().solvent_volume_m3())
            .sum(),
        2.0,
        5e-13,
    );
}

#[test]
fn multicomponent_cross_interactions_preserve_independent_sodium_potassium_exchange() {
    // Independent nested scalar roots eliminate charge exactly:
    // nL=[r,1+t-r,t]. GHK current determines psi; sodium and chloride BE
    // determine r,t; water BE closes x. No production Newton/Jacobian is used
    // to obtain these constants. Equal valence does not imply equal transfer.
    let case = MixtureCase::sodium_potassium_chloride();
    let mut pair = case.build();
    assert_eq!(
        pair.compartments()[0].solution().solute_moles(),
        case.amounts[0].as_slice()
    );
    assert!(pair.endpoint().volume_flux_m_s > 0.0);
    let before = mixture_potential(&pair, &case);
    let report = pair.step(0.2, controls(), || false).unwrap();
    let expected = [
        [
            0.952_776_100_630_362_8,
            0.269_505_027_614_002_97,
            0.222_281_128_244_365_74,
        ],
        [
            0.247_223_899_369_637_18,
            0.530_494_972_385_997_1,
            0.777_718_871_755_634_3,
        ],
    ];
    for (side, cell) in pair.compartments().iter().enumerate() {
        for (i, &n) in cell.solution().solute_moles().iter().enumerate() {
            near(n, expected[side][i], 2e-9);
        }
    }
    let x = report.solvent_transfer_m3;
    near(x, -0.007_439_026_018_840_444_5, 2e-9);
    near(psi(&report.endpoint), -0.311_836_303_951_942_73, 2e-9);
    assert_eq!(report.solute_transfer_mol.len(), 3);
    assert_eq!(report.solute_balance_residual_mol.len(), 3);
    assert_eq!(report.solute_root_residual_mol.len(), 3);
    assert_eq!(report.solute_transfer_mol.into_iter().count(), 3);
    assert!(report.solute_transfer_mol[0] > 0.0);
    assert!(report.solute_transfer_mol[1] < 0.0);
    assert!(report.solute_transfer_mol[2] < 0.0);
    let original_flux = independent_mixture_fluxes(&pair, &case, psi(&report.endpoint));
    for (i, flux) in original_flux.into_iter().enumerate() {
        near(report.solute_transfer_mol[i], 0.2 * flux, 2e-9);
        near(report.solute_root_residual_mol[i], 0.0, 1e-9);
    }
    near(report.endpoint.charge_flux_residual_mol_m2_s, 0.0, 2e-13);
    near(report.water_root_residual_m3, 0.0, 1e-9);
    assert_mixture_conserved(&pair, &case);
    let exact_work = -0.5 * RT * (x * x - (-(x * x)).ln_1p());
    near(report.pressure_work_into_solution_j, exact_work, 3e-9);
    near(
        report.pressure_work_into_solution_j,
        -0.137_185_169_028_683_77,
        2e-6,
    );
    let change = mixture_potential(&pair, &case) - before;
    near(report.total_potential_change_j, change, 3e-8);
    near(change, -494.734_292_876_233_9, 2e-5);
    let dissipation = report.endpoint.total_dissipation_w_m2.unwrap();
    near(dissipation, 2166.492_754_339_36, 2e-4);
    assert!(change <= -0.2 * dissipation + 1e-6);
}

#[test]
fn multicomponent_species_permutation_and_side_reversal_preserve_the_physics() {
    let case = MixtureCase::sodium_potassium_chloride();
    let mut baseline = case.build();
    let original = baseline.step(0.2, controls(), || false).unwrap();
    let order = [2, 0, 1];
    let reordered_case = case.permuted(&order);
    let mut reordered = reordered_case.build();
    let reordered_report = reordered.step(0.2, controls(), || false).unwrap();
    near(
        psi(&reordered_report.endpoint),
        psi(&original.endpoint),
        2e-9,
    );
    near(
        reordered_report.solvent_transfer_m3,
        original.solvent_transfer_m3,
        2e-9,
    );
    for i in 0..3 {
        near(
            reordered_report.solute_transfer_mol[i],
            original.solute_transfer_mol[order[i]],
            2e-9,
        );
        for side in 0..2 {
            near(
                reordered.compartments()[side].solution().solute_moles()[i],
                baseline.compartments()[side].solution().solute_moles()[order[i]],
                2e-9,
            );
        }
    }
    near(
        reordered_report.total_potential_change_j,
        original.total_potential_change_j,
        2e-5,
    );
    assert_mixture_conserved(&reordered, &reordered_case);

    let mut reversed_case = case.clone();
    reversed_case.amounts.swap(0, 1);
    reversed_case.fixed.swap(0, 1);
    reversed_case.trapped.swap(0, 1);
    let mut reversed = reversed_case.build();
    let reverse = reversed.step(0.2, controls(), || false).unwrap();
    near(psi(&reverse.endpoint), -psi(&original.endpoint), 2e-9);
    near(
        reverse.solvent_transfer_m3,
        -original.solvent_transfer_m3,
        2e-9,
    );
    for i in 0..3 {
        near(
            reverse.solute_transfer_mol[i],
            -original.solute_transfer_mol[i],
            2e-9,
        );
        for side in 0..2 {
            near(
                reversed.compartments()[side].solution().solute_moles()[i],
                baseline.compartments()[1 - side].solution().solute_moles()[i],
                2e-9,
            );
        }
    }
    near(
        reverse.pressure_work_into_solution_j,
        original.pressure_work_into_solution_j,
        2e-7,
    );
    near(
        reverse.total_potential_change_j,
        original.total_potential_change_j,
        2e-5,
    );
    assert_mixture_conserved(&reversed, &reversed_case);
}

#[test]
fn full_layout_keeps_neutral_blocked_and_absent_species_distinct() {
    let base = MixtureCase::sodium_potassium_chloride();
    let mut case = base.clone();
    case.amounts[0].extend([0.25, 0.03, 0.0]);
    case.amounts[1].extend([0.05, 0.01, 0.0]);
    case.valences.extend([0, 1, -2]);
    case.permeances.extend([0.25, 0.0, 0.4]);
    case.fixed = [-1.03, -0.01]; // Explicitly neutralize the retained +1 pool.
    for side in 0..2 {
        case.amounts[side].resize(64, 0.0);
    }
    case.valences.resize(64, 0);
    case.permeances.resize(64, 0.0);
    let mut matrix = vec![0.0; 64 * 64];
    for i in 0..3 {
        for j in 0..3 {
            matrix[i * 64 + j] = base.matrix.as_ref().unwrap()[i * 3 + j];
        }
    }
    matrix[3 * 64 + 3] = 0.4;
    matrix[4 * 64 + 4] = 0.25;
    matrix[5 * 64 + 5] = 0.6;
    case.matrix = Some(matrix);
    let mut pair = case.build();
    let before = mixture_potential(&pair, &case);
    let report = pair.step(0.2, controls(), || false).unwrap();
    assert_eq!(pair.species_count(), 64);
    assert_eq!(report.solute_transfer_mol.len(), 64);
    assert_eq!(report.solute_balance_residual_mol.len(), 64);
    assert_eq!(report.solute_root_residual_mol.len(), 64);
    assert_eq!(report.solute_transfer_mol.into_iter().count(), 64);
    near(report.solvent_transfer_m3, -0.025_501_813_021_900_133, 2e-9);
    near(psi(&report.endpoint), -0.278_692_947_477_789_33, 2e-9);
    let expected_left = [
        0.953_515_031_056_831,
        0.270_309_506_093_384_9,
        0.223_824_537_150_215_89,
        0.240_780_913_573_174_66,
        0.03,
        0.0,
    ];
    for (i, &expected) in expected_left.iter().enumerate() {
        near(
            pair.compartments()[0].solution().solute_moles()[i],
            expected,
            2e-9,
        );
    }
    assert!(report.solute_transfer_mol[3] > 0.0); // Neutral urea exchanges.
    for side in 0..2 {
        assert_eq!(
            pair.compartments()[side].solution().solute_moles()[4],
            case.amounts[side][4]
        ); // Blocked charged pool stays exact.
        for i in 5..64 {
            assert_eq!(pair.compartments()[side].solution().solute_moles()[i], 0.0);
        }
    }
    for i in 4..64 {
        assert_eq!(report.solute_transfer_mol[i], 0.0);
    }
    let flux = independent_mixture_fluxes(&pair, &case, psi(&report.endpoint));
    for (i, &j) in flux.iter().enumerate() {
        near(report.solute_transfer_mol[i], 0.2 * j, 2e-9);
    }
    assert_mixture_conserved(&pair, &case);
    near(report.total_potential_change_j, -547.573_113_690_522, 3e-5);
    near(
        report.total_potential_change_j,
        mixture_potential(&pair, &case) - before,
        5e-8,
    );
    near(
        report.pressure_work_into_solution_j,
        -1.612_433_260_740_403,
        2e-6,
    );
    assert!(report.endpoint.total_dissipation_w_m2.unwrap() > 0.0);
}

#[test]
fn all_sixty_four_active_neutral_species_match_independent_analytic_backward_euler() {
    let left: Vec<f64> = (0..64).map(|i| 1.0 + 0.01 * f64::from(i)).collect();
    let permeances: Vec<f64> = (0..64).map(|i| 0.1 + 0.01 * f64::from(i)).collect();
    let case = MixtureCase {
        amounts: [left, vec![0.2; 64]],
        valences: vec![0; 64],
        permeances,
        matrix: None,
        hydraulic: 0.0,
        fixed: [0.0; 2],
        trapped: [0.0; 2],
    };
    let mut pair = case.build();
    let dt = 0.2;
    let report = pair.step(dt, controls(), || false).unwrap();
    assert_eq!(pair.species_count(), 64);
    assert_eq!(report.solute_transfer_mol.len(), 64);
    assert_eq!(report.solute_balance_residual_mol.len(), 64);
    assert_eq!(report.solute_root_residual_mol.len(), 64);
    assert_eq!(report.solute_transfer_mol.into_iter().count(), 64);
    assert_eq!(report.endpoint.potential, DonnanPotential::Indeterminate);
    assert_eq!(report.solvent_transfer_m3, 0.0);
    for i in 0..64 {
        let p = case.permeances[i];
        let expected = dt * p * (case.amounts[0][i] - 0.2) / (1.0 + 2.0 * dt * p);
        near(report.solute_transfer_mol[i], expected, 2e-10);
        near(
            pair.compartments()[0].solution().solute_moles()[i],
            case.amounts[0][i] - expected,
            2e-10,
        );
        near(
            pair.compartments()[1].solution().solute_moles()[i],
            0.2 + expected,
            2e-10,
        );
        near(report.solute_root_residual_mol[i], 0.0, 1e-9);
    }
    assert_mixture_conserved(&pair, &case);
    assert_eq!(report.pressure_work_into_solution_j, 0.0);
    assert!(report.total_potential_change_j < 0.0);
}

#[test]
fn adaptive_accuracy_includes_a_neutral_pool_after_the_charged_species() {
    let mut case = MixtureCase::sodium_potassium_chloride();
    case.amounts[0].push(0.9);
    case.amounts[1].push(0.1);
    case.valences.push(0);
    case.permeances = vec![0.0, 0.0, 0.0, 3.0];
    case.matrix = None;
    case.hydraulic = 0.0;
    case.trapped = [0.0; 2];
    let mut pair = case.build();
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let temporal = IonicAdvanceControls {
        maximum_step_s: 0.5,
        relative_tolerance: 1e-4,
        absolute_amount_tolerance_mol: 1e-10,
        absolute_water_tolerance_m3: 1.0,
        step: controls(),
        ..IonicAdvanceControls::default()
    };
    let report = pair
        .advance(0.5, temporal, &mut workspace, || false)
        .unwrap();
    // Neutral diffusion alone obeys nL=.5+.4 exp(-6t). The water and all
    // charged pools are constant, so omitting the fourth error component
    // would falsely accept one large step and miss this analytic trajectory.
    let exact = 0.5 + 0.4 * (-3.0_f64).exp();
    near(
        pair.compartments()[0].solution().solute_moles()[3],
        exact,
        6e-4,
    );
    assert!(report.accepted_steps > 10);
    assert!(report.maximum_scaled_error <= 1.0);
    assert_eq!(report.solute_transfer_mol.len(), 4);
    assert_eq!(
        report.endpoint.endpoint.potential,
        DonnanPotential::Indeterminate
    );
    for side in 0..2 {
        assert_eq!(
            pair.compartments()[side].solution().solvent_volume_m3(),
            1.0
        );
        for i in 0..3 {
            assert_eq!(
                pair.compartments()[side].solution().solute_moles()[i],
                case.amounts[side][i]
            );
        }
    }
    assert_mixture_conserved(&pair, &case);
    assert_eq!(pair.time_s(), 0.5);
}

#[test]
fn multicomponent_supported_water_seed_keeps_original_concentration_caps() {
    // The third, exactly absent neutral selects the general path while leaving
    // the independent capped binary problem unchanged. At zero water transfer,
    // any finite salt entry exceeds the initially saturated cation cap.
    let caps = [1.2, 1.0, 1.0];
    let model = Arc::new(
        VirialSolutionModel::new(
            T,
            18e-6,
            1.0,
            caps.to_vec(),
            vec![0.8, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.4],
        )
        .unwrap(),
    );
    let amounts = [[1.2, 0.2, 0.0], [0.8, 0.8, 0.0]];
    let cells = [0, 1].map(|side| {
        OsmoticCompartment::new(
            AqueousSolution::nonideal(Arc::clone(&model), 1.0, amounts[side].to_vec()).unwrap(),
            OsmoticPressureVolume::new(1.0, [0.0, RT][side], 0.0).unwrap(),
            0.01,
        )
        .unwrap()
    });
    let membrane = ElectrodiffusiveMembrane::new(
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        vec![1, -1, 0],
        [0.7, 0.2, 0.3]
            .map(|p| SolutePermeance::new(p).unwrap())
            .to_vec(),
    )
    .unwrap();
    let mut pair =
        IonicOsmoticPair::new(cells, membrane, [-1.0, 0.0], [0.0; 2], 1.0, 1e-11, 1 << 20).unwrap();
    let report = pair.step(0.01, controls(), || false).unwrap();
    near(report.solvent_transfer_m3, -0.008_417_969_072_558_327, 2e-9);
    for i in 0..2 {
        near(
            report.solute_transfer_mol[i],
            -0.000_880_660_644_880_319,
            2e-9,
        );
        near(report.solute_root_residual_mol[i], 0.0, 1e-9);
    }
    assert_eq!(report.solute_transfer_mol.len(), 3);
    assert_eq!(report.solute_transfer_mol[2], 0.0);
    near(psi(&report.endpoint), -0.828_210_014_185_601_7, 2e-9);
    for cell in pair.compartments() {
        for (i, &n) in cell.solution().solute_moles().iter().enumerate() {
            assert!(n / cell.solution().solvent_volume_m3() <= caps[i]);
        }
        assert_eq!(cell.solution().solute_moles()[2], 0.0);
    }
    near(
        report.pressure_work_into_solution_j,
        20.867_783_607_393_96,
        2e-6,
    );
    near(
        report.total_potential_change_j,
        -20.503_347_687_013_104,
        2e-6,
    );
}

#[test]
fn multicomponent_late_cancellation_restores_every_pool_and_retries_exactly() {
    let case = MixtureCase::sodium_potassium_chloride();
    let temporal = IonicAdvanceControls {
        maximum_step_s: 0.1,
        relative_tolerance: 5e-3,
        absolute_amount_tolerance_mol: 1e-8,
        absolute_water_tolerance_m3: 1e-8,
        step: controls(),
        ..IonicAdvanceControls::default()
    };
    let mut successful = case.build();
    let mut successful_workspace = successful.adaptive_workspace(1 << 20).unwrap();
    let mut successful_polls = 0usize;
    let accepted = successful
        .advance(0.3, temporal, &mut successful_workspace, || {
            successful_polls += 1;
            false
        })
        .unwrap();
    assert!(accepted.accepted_steps >= 3);
    assert!(accepted.maximum_scaled_error <= 1.0);
    assert!(successful_polls > 10);
    assert_mixture_conserved(&successful, &case);

    let mut pair = case.build();
    let original = snapshot(&pair);
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let mut polls = 0usize;
    let mut cancellation_reached = false;
    let refused = pair.advance(0.3, temporal, &mut workspace, || {
        polls += 1;
        cancellation_reached = polls >= successful_polls - 1;
        cancellation_reached
    });
    assert!(cancellation_reached);
    assert!(matches!(
        refused,
        Err(IonicOsmoticError::Network(OsmoticNetworkError::Cancelled))
    ));
    assert_eq!(snapshot(&pair), original);
    let retry = pair
        .advance(0.3, temporal, &mut workspace, || false)
        .unwrap();
    assert_eq!(snapshot(&pair), snapshot(&successful));
    assert_eq!(retry.accepted_steps, accepted.accepted_steps);
    assert_eq!(retry.rejected_steps, accepted.rejected_steps);
    assert_eq!(retry.solute_transfer_mol, accepted.solute_transfer_mol);
    assert_eq!(
        retry.pressure_work_into_solution_j.to_bits(),
        accepted.pressure_work_into_solution_j.to_bits()
    );
}
