//! G1 analytic limits, G3 conservation/refinement, and atomic refusal regressions.

use super::*;
use fs_material::osmosis::{HydraulicPermeance, ReflectionCoefficient, SolutePermeance};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const TEMPERATURE: f64 = 298.15;
const RT: f64 = R * TEMPERATURE;

fn cell(volume: f64, amounts: &[f64], pressure: f64, elastance: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(TEMPERATURE, 18e-6, volume, amounts.to_vec()).unwrap(),
        OsmoticPressureVolume::new(volume, pressure, elastance).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}

fn membrane(left: usize, right: usize, lp: f64, species: &[(f64, f64)]) -> OsmoticMembrane {
    let law = MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        species
            .iter()
            .map(|&(sigma, permeance)| {
                (
                    ReflectionCoefficient::new(sigma).unwrap(),
                    SolutePermeance::new(permeance).unwrap(),
                )
            })
            .collect(),
    )
    .unwrap();
    OsmoticMembrane::new(left, right, 1.0, law).unwrap()
}

fn network(cells: Vec<OsmoticCompartment>, membranes: Vec<OsmoticMembrane>) -> OsmoticNetwork {
    let species = (0..cells[0].solution().solute_moles().len())
        .map(|i| format!("solute-{i}"))
        .collect();
    OsmoticNetwork::new(species, cells, membranes, 1 << 20).unwrap()
}

fn controls(step: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: step,
        relative_tolerance: 1e-6,
        absolute_volume_tolerance_m3: 1e-10,
        absolute_amount_tolerance_mol: 1e-10,
        maximum_steps: 50_000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.15e}, expected={expected:.15e}, tolerance={tolerance:.3e}"
    );
}

fn snapshot(network: &OsmoticNetwork) -> Vec<u64> {
    let mut bits = vec![network.time_s().to_bits()];
    for cell in network.compartments() {
        bits.push(cell.solution().solvent_volume_m3().to_bits());
        bits.extend(cell.solution().solute_moles().iter().map(|x| x.to_bits()));
    }
    bits
}

#[test]
fn pure_diffusion_into_empty_recipient_matches_unequal_volume_exponential() {
    let mut system = network(
        vec![cell(1.0, &[30.0], 0.0, 0.0), cell(2.0, &[0.0], 0.0, 0.0)],
        vec![membrane(0, 1, 0.0, &[(0.0, 0.2)])],
    );
    let report = system.advance(4.0, controls(0.25), || false).unwrap();
    // Exact two-reservoir diffusion: dc/dt = -AP(1/VL+1/VR) dc.
    let difference = 30.0 * (-0.2_f64 * (1.0 + 0.5) * 4.0).exp();
    let expected_left = 10.0 + (2.0 / 3.0) * difference;
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    near(left.solute_moles()[0], expected_left, 0.01);
    near(right.solute_moles()[0], 30.0 - expected_left, 0.01);
    assert!(right.solute_moles()[0] > 0.0);
    near(
        left.solute_moles()[0] + right.solute_moles()[0],
        30.0,
        2e-12,
    );
    assert_eq!(
        [left.solvent_volume_m3(), right.solvent_volume_m3()],
        [1.0, 2.0]
    );
    assert!(report.maximum_scaled_error <= 1.0);
    assert_eq!(report.water_balance_residual_m3, 0.0);
    near(report.solute_balance_residual_mol[0], 0.0, 2e-12);
}

#[test]
fn compliant_pressure_has_exact_implicit_step_exponential_transient_and_equilibrium() {
    let left = cell(1.0, &[], 200_000.0, 100_000.0);
    let right = cell(1.0, &[], 100_000.0, 100_000.0);
    let link = membrane(0, 1, 1e-6, &[]);
    let endpoint = exchange(
        &left,
        &right,
        &link,
        4.0,
        controls(0.2),
        &mut 0,
        &mut || false,
    )
    .unwrap();
    let implicit_transfer = 4.0 * 1e-6 * 100_000.0 / (1.0 + 4.0 * 1e-6 * 200_000.0);
    near(1.0 - endpoint.volumes[0], implicit_transfer, 1e-10);
    let mut system = network(vec![left, right], vec![link]);
    system.advance(4.0, controls(0.2), || false).unwrap();
    let exact_transfer = 0.5 * (1.0 - (-0.2_f64 * 4.0).exp());
    near(
        1.0 - system.compartments()[0].solution().solvent_volume_m3(),
        exact_transfer,
        4e-4,
    );
    system.advance(80.0, controls(0.2), || false).unwrap();
    near(
        system.compartments()[0].solution().solvent_volume_m3(),
        0.5,
        2e-7,
    );
    near(
        system.compartments()[1].solution().solvent_volume_m3(),
        1.5,
        2e-7,
    );
    near(
        system.compartments()[0].pressure_pa().unwrap(),
        150_000.0,
        0.03,
    );
    near(
        system.compartments()[1].pressure_pa().unwrap(),
        150_000.0,
        0.03,
    );
}

#[test]
fn retained_solutes_reach_one_to_two_volume_ratio_without_creating_water_or_solute() {
    let mut system = network(
        vec![
            cell(1.0, &[50.0], 101_325.0, 0.0),
            cell(1.0, &[100.0], 101_325.0, 0.0),
        ],
        vec![membrane(0, 1, 1e-6, &[(1.0, 0.0)])],
    );
    let report = system.advance(80.0, controls(0.2), || false).unwrap();
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    near(left.solvent_volume_m3(), 2.0 / 3.0, 2e-8);
    near(right.solvent_volume_m3(), 4.0 / 3.0, 2e-8);
    near(
        left.solvent_volume_m3() + right.solvent_volume_m3(),
        2.0,
        2e-12,
    );
    assert_eq!(left.solute_moles(), &[50.0]);
    assert_eq!(right.solute_moles(), &[100.0]);
    near(
        left.osmotic_pressure_pa().unwrap(),
        right.osmotic_pressure_pa().unwrap(),
        0.01,
    );
    near(report.water_balance_residual_m3, 0.0, 2e-12);
    assert_eq!(report.solute_balance_residual_mol[0], 0.0);
}

#[test]
fn nonselective_pure_advection_keeps_uniform_concentration_while_volumes_change() {
    let mut system = network(
        vec![
            cell(1.5, &[60.0], 100_000.0, 0.0),
            cell(0.5, &[20.0], 0.0, 0.0),
        ],
        vec![membrane(0, 1, 1e-6, &[(0.0, 0.0)])],
    );
    system.advance(2.0, controls(0.1), || false).unwrap();
    for (cell, volume) in system.compartments().iter().zip([1.3, 0.7]) {
        near(cell.solution().solvent_volume_m3(), volume, 1e-9);
        near(
            cell.solution().solute_moles()[0] / cell.solution().solvent_volume_m3(),
            40.0,
            1e-8,
        );
    }
}

fn free_energy(system: &OsmoticNetwork) -> f64 {
    system
        .compartments()
        .iter()
        .map(|cell| {
            let volume = cell.solution().solvent_volume_m3();
            let displacement = volume - cell.pressure.reference_volume_m3;
            let mechanical = cell.pressure.reference_pressure_pa * displacement
                + 0.5 * cell.pressure.elastance_pa_per_m3 * displacement * displacement;
            mechanical
                + RT * cell
                    .solution()
                    .solute_moles()
                    .iter()
                    .map(|&amount| {
                        // Independent ideal-solution Helmholtz oracle, including n log n at zero.
                        if amount == 0.0 {
                            0.0
                        } else {
                            amount * ((amount / volume).ln() - 1.0)
                        }
                    })
                    .sum::<f64>()
        })
        .sum()
}

#[test]
fn leaking_mixture_conserves_each_species_and_decreases_isothermal_free_energy() {
    let mut system = network(
        vec![
            cell(1.0, &[20.0, 5.0], 100_000.0, 20_000.0),
            cell(1.0, &[5.0, 30.0], 100_000.0, 20_000.0),
        ],
        vec![membrane(0, 1, 2e-6, &[(0.8, 0.03), (0.2, 0.1)])],
    );
    let mut previous = free_energy(&system);
    for _ in 0..12 {
        let report = system.advance(0.4, controls(0.05), || false).unwrap();
        let current = free_energy(&system);
        assert!(current.is_finite() && current <= previous + 1e-10 * previous.abs());
        previous = current;
        let cells = system.compartments();
        near(
            cells.iter().map(|c| c.solution().solvent_volume_m3()).sum(),
            2.0,
            5e-12,
        );
        for (i, expected) in [25.0, 35.0].into_iter().enumerate() {
            near(
                cells.iter().map(|c| c.solution().solute_moles()[i]).sum(),
                expected,
                5e-12,
            );
            near(report.solute_balance_residual_mol[i], 0.0, 5e-12);
            assert!(cells.iter().all(|c| c.solution().solute_moles()[i] >= 0.0));
        }
    }
    assert!((system.compartments()[0].solution().solvent_volume_m3() - 1.0).abs() > 1e-4);
    assert!((system.compartments()[0].solution().solute_moles()[0] - 20.0).abs() > 1.0);
}

#[test]
fn large_implicit_step_selects_near_origin_flow_direction_not_remote_root() {
    // RT scaling embeds the independently solved dimensionless counterexample
    // in kelvin/Pa units: this residual also has a WRONG positive root ~0.34519.
    let left = cell(1.0, &[100.0], 8.0 * RT, 0.0);
    let right = cell(1.0, &[0.0], 0.0, 0.0);
    let link = membrane(0, 1, 1.0 / RT, &[(0.1, 0.0)]);
    let endpoint = exchange(
        &left,
        &right,
        &link,
        10.0,
        controls(10.0),
        &mut 0,
        &mut || false,
    )
    .unwrap();
    let transfer = 1.0 - endpoint.volumes[0];
    assert!(
        transfer < 0.0,
        "initial hydraulic-minus-osmotic drive is negative"
    );
    near(transfer, -0.24616546375867243, 1e-9);
    near(endpoint.volumes[0] + endpoint.volumes[1], 2.0, 1e-14);
    assert_eq!([endpoint.left[0], endpoint.right[0]], [100.0, 0.0]);
}

fn chain() -> OsmoticNetwork {
    network(
        vec![
            cell(1.0, &[30.0], 0.0, 0.0),
            cell(1.0, &[0.0], 0.0, 0.0),
            cell(1.0, &[0.0], 0.0, 0.0),
        ],
        vec![
            membrane(0, 1, 0.0, &[(0.0, 0.2)]),
            membrane(1, 2, 0.0, &[(0.0, 0.2)]),
        ],
    )
}

#[test]
fn three_cell_graph_conserves_and_refines_to_independent_diffusion_eigenmodes() {
    let t = 3.2;
    let slow = (-0.2_f64 * t).exp();
    let fast = (-0.6_f64 * t).exp();
    let exact = [
        10.0 + 15.0 * slow + 5.0 * fast,
        10.0 - 10.0 * fast,
        10.0 - 15.0 * slow + 5.0 * fast,
    ];
    let mut errors = Vec::new();
    for step in [0.4, 0.2, 0.1] {
        let mut system = chain();
        let mut fixed = controls(step);
        fixed.relative_tolerance = 0.9; // Hold the authored rung rather than adapt it away.
        system.advance(t, fixed, || false).unwrap();
        let amounts: Vec<_> = system
            .compartments()
            .iter()
            .map(|c| c.solution().solute_moles()[0])
            .collect();
        near(amounts.iter().sum(), 30.0, 1e-12);
        assert!(
            system
                .compartments()
                .iter()
                .all(|c| c.solution().solvent_volume_m3() == 1.0)
        );
        errors.push(
            amounts
                .iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "errors={errors:?}"
    );
    assert!(errors[2] < 0.07);
}

#[test]
fn late_cancellation_and_work_exhaustion_preserve_checkpoint_and_allow_exact_retry() {
    let original = chain();
    let before = snapshot(&original);
    let mut config = controls(0.1);
    config.relative_tolerance = 0.9; // First sweeps succeed before the injected interruption.
    let mut cancelled = original.clone();
    let mut polls = 0;
    let error = cancelled
        .advance(1.0, config, || {
            polls += 1;
            polls >= 30
        })
        .unwrap_err();
    assert!(matches!(error, OsmoticNetworkError::Cancelled));
    assert_eq!(snapshot(&cancelled), before);
    let mut exhausted = original.clone();
    let mut bounded = config;
    bounded.maximum_steps = 1;
    assert!(matches!(
        exhausted.advance(1.0, bounded, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(snapshot(&exhausted), before);
    let mut reference = original;
    reference.advance(1.0, config, || false).unwrap();
    cancelled.advance(1.0, config, || false).unwrap();
    exhausted.advance(1.0, config, || false).unwrap();
    assert_eq!(snapshot(&cancelled), snapshot(&reference));
    assert_eq!(snapshot(&exhausted), snapshot(&reference));
}

#[test]
fn almost_complete_advection_preserves_finite_donor_concentration_in_both_directions() {
    for water in [1.0_f64.next_down(), -1.0_f64.next_down()] {
        let [left, right] = [water, 0.0].map(|pressure| {
            OsmoticCompartment::new(
                IdealDiluteSolution::new(TEMPERATURE, 18e-6, 1.0, vec![1.0]).unwrap(),
                OsmoticPressureVolume::new(1.0, pressure, 0.0).unwrap(),
                1e-30,
            )
            .unwrap()
        });
        let link = membrane(0, 1, 1.0, &[(0.0, 0.0)]);
        let endpoint = pair_endpoint(&left, &right, &link, 1.0, water).unwrap();
        // Subtracting a rounded, nearly complete transfer used to erase the
        // donor's representable remainder and its order-one concentration.
        assert!(endpoint.left[0] > 0.0 && endpoint.right[0] > 0.0);
        near(endpoint.left[0] / endpoint.volumes[0], 1.0, 1e-15);
        near(endpoint.right[0] / endpoint.volumes[1], 1.0, 1e-15);
        assert_eq!(endpoint.left[0] + endpoint.right[0], 2.0);
        near(endpoint.residual_m3, 0.0, 1e-15);
    }
}

#[test]
fn inadmissible_distant_volume_floor_does_not_block_a_valid_nearby_root() {
    // Extreme admitted SI inputs exercise arithmetic admission, not a claim
    // that the ideal-dilute model describes concentrated physical solutions.
    let [left, right] = [1e305, 0.0].map(|pressure| {
        OsmoticCompartment::new(
            IdealDiluteSolution::new(TEMPERATURE, 18e-6, 1.0, vec![1e300]).unwrap(),
            OsmoticPressureVolume::new(1.0, pressure, 1e307).unwrap(),
            1e-300,
        )
        .unwrap()
    });
    let link = membrane(0, 1, 1e-307, &[(1.0, 0.0)]);
    let old_far_bound = (left.solution().solvent_volume_m3() - left.minimum_volume_m3).next_down();
    assert!(pair_endpoint(&left, &right, &link, 1.0, old_far_bound).is_err());
    let endpoint = exchange(
        &left,
        &right,
        &link,
        1.0,
        controls(1.0),
        &mut 0,
        &mut || false,
    )
    .unwrap();
    let water = 1.0 - endpoint.volumes[0];
    assert!(water > 0.003 && water < 0.004, "nearby transfer={water}");
    assert_eq!([endpoint.left[0], endpoint.right[0]], [1e300, 1e300]);
    near(endpoint.volumes[0] + endpoint.volumes[1], 2.0, 1e-14);
    near(endpoint.residual_m3, 0.0, 2e-12);
}

#[test]
fn individually_finite_compartments_with_overflowing_network_water_total_refuse() {
    let cells = vec![cell(1e308, &[], 0.0, 0.0), cell(1e308, &[], 0.0, 0.0)];
    assert!(matches!(
        OsmoticNetwork::new(vec![], cells, vec![], 1 << 20),
        Err(OsmoticNetworkError::Invalid(
            "network inventory total is not representable"
        ))
    ));
}

#[test]
fn empty_edge_network_observes_temporal_cancellation_before_step_budget_exhaustion() {
    let mut system = network(vec![cell(1.0, &[], 0.0, 0.0)], vec![]);
    let before = snapshot(&system);
    let mut config = controls(1.0);
    config.maximum_steps = 1;
    let mut polls = 0;
    let result = system.advance(2.0, config, || {
        polls += 1;
        polls == 2
    });
    assert!(matches!(result, Err(OsmoticNetworkError::Cancelled)));
    assert_eq!(polls, 2);
    assert_eq!(snapshot(&system), before);
}
