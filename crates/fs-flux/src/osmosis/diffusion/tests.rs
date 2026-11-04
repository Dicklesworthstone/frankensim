use super::super::{OsmoticPressureVolume, OsmoticStepControls};
use super::*;
use fs_material::osmosis::{
    AqueousSolution,
    nonideal::{NonidealSolutionError, VirialSolutionModel},
};
use fs_qty::MolarMass;
use fs_thermochem::{SpeciesId, transport::maxwell_stefan::MaxwellStefanMolarSolveV1};
use std::sync::Arc;

mod electrical;
mod implicit;
mod multicomponent;
mod source;

const T: f64 = 298.15;
const VBAR: f64 = 18e-6;
const VOLUME: f64 = 1e-6;

fn source(reference: f64, interactions: Vec<f64>, caps: Vec<f64>) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, VBAR, reference, caps, interactions).unwrap())
}
fn cell(model: &Arc<VirialSolutionModel>, concentrations: &[f64]) -> OsmoticCompartment {
    let solution = AqueousSolution::nonideal(
        model.clone(),
        VOLUME,
        concentrations.iter().map(|c| c * VOLUME).collect(),
    )
    .unwrap();
    OsmoticCompartment::new(
        solution,
        OsmoticPressureVolume::new(VOLUME, 100000.0, 1e8).unwrap(),
        0.1 * VOLUME,
    )
    .unwrap()
}
fn physical() -> SolutionDiffusionControls {
    SolutionDiffusionControls {
        constitutive: SolutionDiffusionSolve {
            friction: MaxwellStefanMolarSolveV1::default(),
            entropy_tolerance_w_m3_k: 1e-9,
        },
        maximum_species_fraction: 0.2,
        amount_balance_tolerance_mol: 1e-16,
        free_energy_tolerance_j: 1e-12,
    }
}
fn numerical(relative: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: 50.0,
        relative_tolerance: relative,
        absolute_amount_tolerance_mol: 1e-14,
        absolute_volume_tolerance_m3: 1e-18,
        maximum_steps: 10000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}
fn law(model: &Arc<VirialSolutionModel>) -> SolutionMaxwellStefan {
    let ids = (0..model.species_count())
        .map(|i| {
            (
                SpeciesId::new(format!("solute-{i}")).unwrap(),
                MolarMass::new(0.04 + 0.02 * i as f64),
            )
        })
        .collect();
    let matrix = if model.species_count() == 1 {
        vec![0.0, 1e-9, 1e-9, 0.0]
    } else {
        vec![0.0, 1e-9, 3e-9, 1e-9, 0.0, 5e-10, 3e-9, 5e-10, 0.0]
    };
    SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        ids,
        matrix,
    )
    .unwrap()
}
fn network(
    model: &Arc<VirialSolutionModel>,
    concentrations: &[&[f64]],
    faces: &[(usize, usize, f64)],
) -> OsmoticNetwork {
    let names = (0..model.species_count())
        .map(|i| format!("solute-{i}"))
        .collect();
    let cells = concentrations.iter().map(|c| cell(model, c)).collect();
    let links = faces
        .iter()
        .map(|&(a, b, area)| SolutionDiffusionLink::new(a, b, area, 0.001, law(model)).unwrap())
        .collect();
    OsmoticNetwork::with_solution_diffusion(names, cells, links, physical(), 1 << 22).unwrap()
}
fn concentration(system: &OsmoticNetwork, cell: usize, species: usize) -> f64 {
    system.compartments()[cell].solution().solute_moles()[species]
        / system.compartments()[cell].solution().solvent_volume_m3()
}
fn bits(system: &OsmoticNetwork) -> Vec<u64> {
    let mut values = vec![system.time_s().to_bits()];
    for cell in system.compartments() {
        values.push(cell.solution().solvent_volume_m3().to_bits());
        values.extend(cell.solution().solute_moles().iter().map(|x| x.to_bits()));
    }
    values
}

#[test]
fn g1_binary_window_refines_to_independent_nonlinear_solvent_frame_relaxation() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    // Independent scalar ODE in the stationary-solvent frame. The logarithmic
    // chain rule gives N=D*(cL-cR)*(1+c_log*Vbar)/L, retaining the finite-solvent
    // correction rather than treating a Maxwell--Stefan coefficient as Fick D.
    let rate = |left: f64| {
        let right = 120.0 - left;
        let mean = (left - right) / (left / right).ln();
        -0.01 / VOLUME * 1e-9 * (left - right) * (1.0 + mean * VBAR) / 0.001
    };
    let mut exact = 100.0;
    let h = 50.0 / 10000.0;
    for _ in 0..10000 {
        let a = rate(exact);
        let b = rate(exact + 0.5 * h * a);
        let c = rate(exact + 0.5 * h * b);
        let d = rate(exact + h * c);
        exact += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for relative in [1e-3, 6.25e-5] {
        let mut system = network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)]);
        let report = system.advance(50.0, numerical(relative), || false).unwrap();
        assert_eq!(system.time_s(), 50.0);
        assert!(report.rejected_steps > 0 && report.accepted_steps > 1);
        assert!(report.maximum_scaled_error <= 1.0);
        assert!(report.residual_evaluations > 3 * report.accepted_steps);
        assert_eq!(report.water_balance_residual_m3, 0.0);
        assert!(
            report.solute_balance_residual_mol[0].abs() <= physical().amount_balance_tolerance_mol
        );
        assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
        assert!(
            (concentration(&system, 0, 0) + concentration(&system, 1, 0) - 120.0).abs() < 1e-10
        );
        for cell in system.compartments() {
            assert_eq!(cell.solution().solvent_volume_m3(), VOLUME);
        }
        errors.push((concentration(&system, 0, 0) - exact).abs());
    }
    assert!(errors[0] < 0.5 && errors[1] < 0.4 * errors[0], "{errors:?}");
}

#[test]
fn g1_g3_actual_nonideal_inventory_moves_uphill_and_lowers_mixing_free_energy() {
    let model = source(1.0, vec![0.001, 0.008, 0.008, 0.1], vec![200.0; 2]);
    let mut system = network(&model, &[&[20.0, 100.0], &[21.0, 1.0]], &[(0, 1, 1e-4)]);
    let pressure = system.compartments()[1]
        .solution()
        .osmotic_pressure_pa()
        .unwrap();
    let report = system.advance(10.0, numerical(1e-4), || false).unwrap();
    assert!(concentration(&system, 0, 0) < 20.0);
    assert!(concentration(&system, 1, 0) > 21.0);
    assert!(concentration(&system, 0, 1) < 100.0);
    assert_ne!(
        system.compartments()[1]
            .solution()
            .osmotic_pressure_pa()
            .unwrap(),
        pressure
    );
    assert!(report.diffusion_free_energy_change_j.unwrap() < -1e-4);
    for residual in &report.solute_balance_residual_mol[..2] {
        assert!(residual.abs() <= physical().amount_balance_tolerance_mol);
    }
}

#[test]
fn g3_branched_face_order_and_thermodynamic_reference_do_not_change_the_window() {
    let ordinary = source(1.0, vec![0.0], vec![200.0]);
    let shifted = source(1e200, vec![0.0], vec![200.0]);
    let mut a = network(
        &ordinary,
        &[&[100.0], &[20.0], &[20.0]],
        &[(0, 1, 0.01), (0, 2, 0.01)],
    );
    let mut b = network(
        &shifted,
        &[&[100.0], &[20.0], &[20.0]],
        &[(2, 0, 0.01), (1, 0, 0.01)],
    );
    let x = a.advance(20.0, numerical(1e-4), || false).unwrap();
    let y = b.advance(20.0, numerical(1e-4), || false).unwrap();
    assert_eq!(bits(&a), bits(&b));
    assert_eq!(x.accepted_steps, y.accepted_steps);
    assert_eq!(x.residual_evaluations, y.residual_evaluations);
    assert_eq!(
        x.diffusion_free_energy_change_j,
        y.diffusion_free_energy_change_j
    );
    assert!(concentration(&a, 0, 0) < 100.0);
    assert!(concentration(&a, 1, 0) > 20.0 && concentration(&a, 2, 0) > 20.0);
}

#[test]
fn g0_g4_source_boundaries_and_late_failures_publish_no_state_or_clock() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut system = network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)]);
    let saved = bits(&system);
    assert!(matches!(
        system.advance(1e-20, numerical(1e-4), || false),
        Err(OsmoticNetworkError::StepUnderflow)
    ));
    assert_eq!(bits(&system), saved);
    let limited = OsmoticStepControls {
        maximum_step_s: 0.5,
        maximum_steps: 1,
        ..numerical(1e-3)
    };
    assert!(matches!(
        system.advance(10.0, limited, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&system), saved);
    let mut probe = system.clone();
    let mut polls = 0;
    probe
        .advance(10.0, numerical(1e-4), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut called = 0;
    assert!(matches!(
        system.advance(10.0, numerical(1e-4), || {
            called += 1;
            called == polls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&system), saved);
    // Actual uphill transport immediately leaves the authored concentration
    // domain. Smaller-step retries must not turn this physical refusal into a
    // work-budget error or silently extend the source beyond its stated cap.
    let bounded = source(1.0, vec![0.001, 0.008, 0.008, 0.1], vec![21.0, 200.0]);
    let mut outside = network(&bounded, &[&[20.0, 100.0], &[21.0, 1.0]], &[(0, 1, 1e-4)]);
    let saved = bits(&outside);
    assert!(matches!(
        outside.advance(0.1, numerical(1e-4), || false),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::Thermodynamics(OsmosisError::Nonideal(
                NonidealSolutionError::OutsideDomain { species: 0, .. }
            ))
        ))
    ));
    assert_eq!(bits(&outside), saved);
}

#[test]
fn g0_closed_face_finishes_time_and_preparation_refuses_unsupported_layouts() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut closed = network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.0)]);
    let saved = bits(&closed);
    let report = closed
        .advance(
            1.25,
            OsmoticStepControls {
                maximum_step_s: 0.5,
                ..numerical(1e-4)
            },
            || false,
        )
        .unwrap();
    assert_eq!(closed.time_s(), 1.25);
    assert_eq!(&bits(&closed)[1..], &saved[1..]);
    assert_eq!(report.residual_evaluations, 0);
    assert_eq!(report.diffusion_free_energy_change_j, Some(0.0));
    let make = |names: Vec<String>, zero: bool, duplicate: bool, budget| {
        let cells = vec![
            cell(&model, &[100.0]),
            cell(&model, &[if zero { 0.0 } else { 20.0 }]),
        ];
        let mut links = vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()];
        if duplicate {
            links.push(links[0].clone());
        }
        OsmoticNetwork::with_solution_diffusion(names, cells, links, physical(), budget)
    };
    assert!(make(vec!["foreign".into()], false, false, 1 << 22).is_err());
    assert!(make(vec!["solute-0".into()], true, false, 1 << 22).is_ok());
    assert!(make(vec!["solute-0".into()], false, true, 1 << 22).is_err());
    assert!(matches!(
        make(vec!["solute-0".into()], false, false, 1),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
}

#[test]
fn g1_pure_solvent_uptake_refines_to_dilute_analytical_diffusion() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let initial = 1e-3;
    let duration = 50.0;
    // Two equal, well-mixed volumes coupled by A*D/L. The dilute-limit
    // concentration difference decays as exp(-2*A*D*t/(V*L)). The retained
    // finite-solvent correction is bounded here by Vbar*initial < 2e-8.
    let rate = 2.0 * 0.01 * 1e-9 / (VOLUME * 0.001);
    let expected = 0.5 * initial * (1.0 - (-rate * duration).exp());
    let mut errors = Vec::new();
    for relative in [1e-3, 1e-5] {
        let mut system = network(&model, &[&[initial], &[0.0]], &[(0, 1, 0.01)]);
        assert_eq!(concentration(&system, 1, 0), 0.0);
        let report = system
            .advance(
                duration,
                OsmoticStepControls {
                    absolute_amount_tolerance_mol: 1e-20,
                    ..numerical(relative)
                },
                || false,
            )
            .unwrap();
        let received = concentration(&system, 1, 0);
        assert!(received > 0.0 && received < 0.5 * initial);
        assert!((concentration(&system, 0, 0) + received - initial).abs() < 1e-16);
        assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
        assert!(report.solute_balance_residual_mol[0].abs() < 1e-20);
        assert_eq!(report.water_balance_residual_m3, 0.0);
        assert_eq!(system.time_s(), duration);
        errors.push((received - expected).abs());
    }
    assert!(errors[0] < initial * 0.01, "{errors:?}");
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
}

#[test]
fn g0_g3_empty_species_stays_absent_while_a_front_crosses_a_network() {
    let build = |reference, reversed| {
        let model = source(reference, vec![0.001, 0.008, 0.008, 0.1], vec![200.0; 2]);
        let faces = if reversed {
            vec![(2, 1, 0.01), (1, 0, 0.01)]
        } else {
            vec![(0, 1, 0.01), (1, 2, 0.01)]
        };
        network(&model, &[&[100.0, 0.0], &[0.0, 0.0], &[0.0, 0.0]], &faces)
    };
    let mut system = build(1.0, false);
    let mut reordered = build(1e200, true);
    let report = system.advance(20.0, numerical(1e-4), || false).unwrap();
    let other = reordered.advance(20.0, numerical(1e-4), || false).unwrap();
    assert_eq!(bits(&system), bits(&reordered));
    assert_eq!(
        report.diffusion_free_energy_change_j,
        other.diffusion_free_energy_change_j
    );
    assert!(concentration(&system, 0, 0) < 100.0);
    assert!(concentration(&system, 1, 0) > 0.0);
    assert!(concentration(&system, 2, 0) > 0.0);
    let total: f64 = (0..3).map(|i| concentration(&system, i, 0)).sum();
    assert!((total - 100.0).abs() < 1e-10);
    for i in 0..3 {
        assert_eq!(concentration(&system, i, 1), 0.0);
        assert_eq!(
            system.compartments()[i].solution().solvent_volume_m3(),
            VOLUME
        );
    }
    assert_eq!(report.solute_balance_residual_mol[1], 0.0);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);

    // With no solute anywhere, every face is still evaluated but the complete
    // physical inventory stays exactly empty; no 0/0 withdrawal is manufactured.
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut empty = network(&model, &[&[0.0], &[0.0]], &[(0, 1, 0.01)]);
    let report = empty.advance(1.0, numerical(1e-4), || false).unwrap();
    assert_eq!(concentration(&empty, 0, 0), 0.0);
    assert_eq!(concentration(&empty, 1, 0), 0.0);
    assert_eq!(report.diffusion_free_energy_change_j, Some(0.0));
}

#[test]
fn g4_zero_inventory_window_rolls_back_late_failure_and_retries_cleanly() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut system = network(&model, &[&[100.0], &[0.0]], &[(0, 1, 0.01)]);
    let original = bits(&system);
    let controls = numerical(1e-4);
    let mut accepted = system.clone();
    let mut polls = 0;
    accepted
        .advance(20.0, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        system.advance(20.0, controls, || {
            visited += 1;
            visited == polls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&system), original);
    assert!(matches!(
        system.advance(
            20.0,
            OsmoticStepControls {
                maximum_steps: 1,
                maximum_step_s: 0.01,
                ..controls
            },
            || false,
        ),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&system), original);
    system.advance(20.0, controls, || false).unwrap();
    assert_eq!(bits(&system), bits(&accepted));
}
