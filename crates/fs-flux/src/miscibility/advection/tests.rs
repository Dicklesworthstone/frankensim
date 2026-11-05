//! G0/G2 prescribed-flow transport: analytic donors, wave refinement and inventories.

use super::{BinaryMiscibilityAdvectionControls, BinaryMiscibilityAdvectionFace};
use crate::miscibility::{
    BinaryMiscibilityCell, BinaryMiscibilityError, BinaryMiscibilityFace, BinaryMiscibilityNetwork,
    BinaryMiscibilityState,
};
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn network(
    chi: f64,
    scaled_kappa: f64,
    concentration: f64,
    volumes: &[f64],
    connections: &[[usize; 2]],
) -> BinaryMiscibilityNetwork {
    let pair = BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        chi,
        0.0,
        Temperature::new(250.0),
        Temperature::new(600.0),
    )
    .unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(concentration),
        scaled_kappa * R * 300.0,
    )
    .unwrap();
    BinaryMiscibilityNetwork::new(
        model,
        volumes
            .iter()
            .map(|&volume_m3| BinaryMiscibilityCell { volume_m3 })
            .collect(),
        connections
            .iter()
            .map(|&cells| BinaryMiscibilityFace {
                cells,
                area_m2: 1.0,
                distance_m: 0.5 * (volumes[cells[0]] + volumes[cells[1]]),
                diffusivity_m2_s: 0.0,
            })
            .collect(),
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
}

fn cycle(volumes: &[f64], concentration: f64) -> BinaryMiscibilityNetwork {
    let connections = (0..volumes.len())
        .map(|i| [i, (i + 1) % volumes.len()])
        .collect::<Vec<_>>();
    network(0.0, 0.0, concentration, volumes, &connections)
}

fn cycle_rates(net: &BinaryMiscibilityNetwork, rate: f64) -> Vec<BinaryMiscibilityAdvectionFace> {
    let wrap = [0, net.cells().len() - 1];
    net.faces()
        .map(|face| {
            let reverse = face.cells == wrap;
            BinaryMiscibilityAdvectionFace {
                forward_m3_s: if reverse { 0.0 } else { rate },
                reverse_m3_s: if reverse { rate } else { 0.0 },
            }
        })
        .collect()
}

fn state(net: &BinaryMiscibilityNetwork, fractions: &[f64]) -> BinaryMiscibilityState {
    net.state_from_mole_fractions(fractions.iter().map(|&x| [x, 1.0 - x]).collect(), || false)
        .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}

#[test]
fn counterflow_transports_pure_donors_even_when_net_volume_flux_is_zero() {
    let net = network(0.0, 0.0, 1.0, &[1.0, 1.0], &[[0, 1]]);
    let initial = state(&net, &[1.0, 0.0]);
    let original = initial.clone();
    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 0.4,
                reverse_m3_s: 0.4,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let candidate = net
        .advect_candidate(
            &initial,
            &flow,
            0.5,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    for (actual, expected) in candidate
        .state()
        .component_moles()
        .iter()
        .zip([[0.8, 0.2], [0.2, 0.8]])
    {
        for species in 0..2 {
            near(actual[species], expected[species], 2e-15);
        }
    }
    assert_eq!(candidate.step().substeps, 1);
    near(candidate.step().maximum_courant, 0.2, 2e-15);
    near(
        candidate.step().transferred_component_moles[0][0],
        0.2,
        2e-15,
    );
    near(
        candidate.step().transferred_component_moles[0][1],
        -0.2,
        2e-15,
    );
    assert_eq!(candidate.state().elapsed_s(), 0.5);
    assert_eq!(initial, original);
}

#[test]
fn a_closed_cycle_uses_the_same_old_donor_state_on_every_face() {
    let net = cycle(&[1.0; 3], 1.0);
    let initial = state(&net, &[1.0, 0.0, 0.25]);
    let flow = net
        .prepare_advection(cycle_rates(&net, 0.5), 0.0, || false)
        .unwrap();
    let candidate = net
        .advect_candidate(
            &initial,
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    // One simultaneous donor step: x_i' = (x_i + x_{i-1}) / 2.
    for (actual, expected) in candidate
        .state()
        .mole_fractions()
        .iter()
        .zip([0.625, 0.5, 0.125])
    {
        near(actual[0], expected, 2e-15);
        near(actual[1], 1.0 - expected, 2e-15);
    }
    assert_eq!(candidate.step().substeps, 1);
}

#[test]
fn periodic_translation_converges_to_the_prescribed_velocity_at_fixed_courant() {
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let h = 1.0 / n as f64;
        let net = cycle(&vec![h; n], 1.0);
        let initial = state(
            &net,
            &(0..n)
                .map(|i| 0.5 + 0.2 * (core::f64::consts::TAU * (i as f64 + 0.5) * h).sin())
                .collect::<Vec<_>>(),
        );
        let flow = net
            .prepare_advection(cycle_rates(&net, 1.0), 0.0, || false)
            .unwrap();
        let candidate = net
            .advect_candidate(
                &initial,
                &flow,
                0.25,
                BinaryMiscibilityAdvectionControls::default(),
                || false,
            )
            .unwrap();
        // Exact periodic transport x(s,t)=1/2+0.2 sin(2 pi (s-t)), u=1 m/s.
        let error = candidate
            .state()
            .mole_fractions()
            .iter()
            .enumerate()
            .map(|(i, x)| {
                let exact =
                    0.5 + 0.2 * (core::f64::consts::TAU * ((i as f64 + 0.5) * h - 0.25)).sin();
                h * (x[0] - exact).abs()
            })
            .sum::<f64>();
        errors.push(error);
        near(candidate.step().maximum_courant, 0.5, 2e-15);
        near(candidate.state().elapsed_s(), 0.25, 2e-15);
    }
    assert!(errors[1] < errors[0] / 1.8, "{errors:?}");
    assert!(errors[2] < errors[1] / 1.8, "{errors:?}");
    assert!(errors[2] < 0.006, "{errors:?}");
}

#[test]
fn unequal_volumes_conserve_moles_and_every_integrated_face_receipt_closes() {
    let net = cycle(&[1.0, 2.0, 4.0], 2.0);
    let initial = state(&net, &[0.2, 0.5, 0.8]);
    let flow = net
        .prepare_advection(cycle_rates(&net, 0.75), 0.0, || false)
        .unwrap();
    let candidate = net
        .advect_candidate(
            &initial,
            &flow,
            4.0,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.step().substeps, 6);
    let mut reconstructed = initial.component_moles().to_vec();
    for (face, transfer) in net
        .faces()
        .zip(&candidate.step().transferred_component_moles)
    {
        let [i, j] = face.cells;
        for species in 0..2 {
            reconstructed[i][species] -= transfer[species];
            reconstructed[j][species] += transfer[species];
        }
    }
    for (i, actual) in candidate.state().component_moles().iter().enumerate() {
        for species in 0..2 {
            near(actual[species], reconstructed[i][species], 3e-14);
        }
        near(actual[0] + actual[1], net.capacities_mol()[i], 3e-14);
    }
    for species in 0..2 {
        near(
            candidate
                .state()
                .component_moles()
                .iter()
                .map(|n| n[species])
                .sum(),
            initial.component_moles().iter().map(|n| n[species]).sum(),
            3e-14,
        );
        near(candidate.step().component_residual_mol[species], 0.0, 3e-14);
    }
    near(candidate.step().maximum_capacity_residual_mol, 0.0, 3e-14);
    // An unweighted concentration sum is not the conserved quantity here.
    let unweighted: f64 = candidate
        .state()
        .mole_fractions()
        .iter()
        .map(|x| x[0])
        .sum();
    assert!((unweighted - 1.5).abs() > 0.01);
}

#[test]
fn divergence_free_flow_keeps_homogeneous_and_absent_species_homogeneous() {
    let net = cycle(&[0.5, 1.0, 2.0], 1.0);
    let flow = net
        .prepare_advection(cycle_rates(&net, 0.25), 0.0, || false)
        .unwrap();
    for fraction in [0.25, 1.0, 0.0] {
        let initial = state(&net, &[fraction; 3]);
        let candidate = net
            .advect_candidate(
                &initial,
                &flow,
                7.0,
                BinaryMiscibilityAdvectionControls::default(),
                || false,
            )
            .unwrap();
        for x in candidate.state().mole_fractions() {
            near(x[0], fraction, 2e-14);
            near(x[1], 1.0 - fraction, 2e-14);
            if fraction == 0.0 || fraction == 1.0 {
                assert_eq!(*x, [fraction, 1.0 - fraction]);
            }
        }
    }
}

#[test]
fn real_trace_species_survive_when_the_dominant_fraction_rounds_to_one() {
    let net = cycle(&[1.0; 3], 1.0);
    let initial = net
        .state_from_mole_fractions(vec![[1.0, 1e-30], [1.0, 3e-30], [1.0, 2e-30]], || false)
        .unwrap();
    let flow = net
        .prepare_advection(cycle_rates(&net, 0.25), 0.0, || false)
        .unwrap();
    let candidate = net
        .advect_candidate(
            &initial,
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    for (actual, expected) in candidate
        .state()
        .component_moles()
        .iter()
        .zip([1.25e-30, 2.5e-30, 2.25e-30])
    {
        assert!(actual[1] > 0.0);
        near(actual[1], expected, expected * 2e-14);
    }
    near(
        candidate
            .state()
            .component_moles()
            .iter()
            .map(|n| n[1])
            .sum(),
        6e-30,
        2e-43,
    );
}

#[test]
fn prescribed_mixing_can_raise_nonconvex_free_energy_without_a_false_refusal() {
    let net = network(4.0, 1e-4, 1.0, &[1.0; 2], &[[0, 1]]);
    let initial = state(&net, &[0.02, 0.98]);
    let flow = net
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 0.2,
                reverse_m3_s: 0.2,
            }],
            0.0,
            || false,
        )
        .unwrap();
    let candidate = net
        .advect_candidate(
            &initial,
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            || false,
        )
        .unwrap();
    let report = candidate.step();
    assert!(report.free_energy_change_j > 100.0);
    near(
        report.free_energy_after_j - report.free_energy_before_j,
        report.free_energy_change_j,
        1e-9,
    );
    near(
        net.free_energy_j(candidate.state(), || false).unwrap(),
        report.free_energy_after_j,
        1e-10,
    );
}

#[test]
fn field_and_flow_require_the_original_graph_and_thermodynamic_owner() {
    let net = cycle(&[1.0; 3], 1.0);
    let other = cycle(&[1.0, 2.0, 1.0], 1.0);
    let other_model = network(0.5, 0.0, 1.0, &[1.0; 3], &[[0, 1], [1, 2], [2, 0]]);
    let initial = state(&net, &[0.2, 0.5, 0.8]);
    let other_state = state(&other, &[0.2, 0.5, 0.8]);
    let flow = net
        .prepare_advection(cycle_rates(&net, 0.25), 0.0, || false)
        .unwrap();
    let other_flow = other
        .prepare_advection(cycle_rates(&other, 0.25), 0.0, || false)
        .unwrap();
    let other_model_flow = other_model
        .prepare_advection(cycle_rates(&other_model, 0.25), 0.0, || false)
        .unwrap();
    for (source, prepared) in [
        (&initial, &other_flow),
        (&other_state, &flow),
        (&initial, &other_model_flow),
    ] {
        assert!(matches!(
            net.advect_candidate(
                source,
                prepared,
                1.0,
                BinaryMiscibilityAdvectionControls::default(),
                || false,
            ),
            Err(BinaryMiscibilityError::Invalid(_))
        ));
    }
}

#[test]
fn preparation_refuses_volume_sources_and_invalid_directional_rates() {
    let net = network(0.0, 0.0, 1.0, &[1.0; 2], &[[0, 1]]);
    assert!(matches!(
        net.prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: 0.5,
                reverse_m3_s: 0.0,
            }],
            1e-12,
            || false,
        ),
        Err(BinaryMiscibilityError::VolumeDivergence { cell: 0, .. })
    ));
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            net.prepare_advection(
                vec![BinaryMiscibilityAdvectionFace {
                    forward_m3_s: invalid,
                    reverse_m3_s: 0.5,
                }],
                0.0,
                || false,
            ),
            Err(BinaryMiscibilityError::Invalid(_))
        ));
    }
    assert!(matches!(
        net.prepare_advection(Vec::new(), 0.0, || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
}

#[test]
fn cancellation_and_substep_exhaustion_publish_no_partial_transport() {
    let net = cycle(&[1.0; 3], 1.0);
    let initial = state(&net, &[0.2, 0.5, 0.8]);
    let original = initial.clone();
    let flow = net
        .prepare_advection(cycle_rates(&net, 1.0), 0.0, || false)
        .unwrap();
    assert!(matches!(
        net.advect_candidate(
            &initial,
            &flow,
            2.0,
            BinaryMiscibilityAdvectionControls {
                maximum_substeps: 2,
                ..BinaryMiscibilityAdvectionControls::default()
            },
            || false,
        ),
        Err(BinaryMiscibilityError::AdvectionWorkLimit {
            required_substeps: 4,
            maximum_substeps: 2,
        })
    ));
    let mut polls = 0;
    assert!(matches!(
        net.advect_candidate(
            &initial,
            &flow,
            2.0,
            BinaryMiscibilityAdvectionControls::default(),
            || {
                polls += 1;
                polls > 16
            },
        ),
        Err(BinaryMiscibilityError::Cancelled)
    ));
    assert!(polls > 16);
    assert_eq!(initial, original);
}
