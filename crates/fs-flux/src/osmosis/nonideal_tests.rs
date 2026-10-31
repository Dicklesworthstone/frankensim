//! Nonideal finite-water transport, independent thermodynamic references, and refusals.
//! Coefficients are authored analytical fixtures, not measured solution calibrations.

use super::*;
use fs_material::osmosis::nonideal::VirialSolutionModel;
use fs_material::osmosis::polarization::FilmPermeance;
use fs_material::osmosis::{HydraulicPermeance, ReflectionCoefficient, SolutePermeance};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::sync::Arc;

const T: f64 = 298.15;
const RT: f64 = R * T;

fn model(caps: Vec<f64>, interaction: Vec<f64>) -> Arc<VirialSolutionModel> {
    Arc::new(VirialSolutionModel::new(T, 18e-6, 1.0, caps, interaction).unwrap())
}

fn ideal(volume: f64, amounts: &[f64]) -> AqueousSolution {
    IdealDiluteSolution::new(T, 18e-6, volume, amounts.to_vec())
        .unwrap()
        .into()
}

fn cell(solution: impl Into<AqueousSolution>, pressure: f64, elastance: f64) -> OsmoticCompartment {
    let solution = solution.into();
    let volume = solution.solvent_volume_m3();
    OsmoticCompartment::new(
        solution,
        OsmoticPressureVolume::new(volume, pressure, elastance).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}

fn law(species: &[(f64, f64)]) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(1.0 / RT).unwrap(),
        species
            .iter()
            .map(|&(sigma, p)| {
                (
                    ReflectionCoefficient::new(sigma).unwrap(),
                    SolutePermeance::new(p).unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
}

fn retained(count: usize) -> OsmoticMembrane {
    OsmoticMembrane::new(0, 1, 1.0, law(&vec![(1.0, 0.0); count])).unwrap()
}

fn names(count: usize) -> Vec<String> {
    (0..count).map(|i| format!("particle-{i}")).collect()
}

fn network(cells: Vec<OsmoticCompartment>) -> OsmoticNetwork {
    let count = cells[0].solution().solute_moles().len();
    OsmoticNetwork::new(names(count), cells, vec![retained(count)], 1 << 20).unwrap()
}

fn controls(step: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: step,
        relative_tolerance: 1e-6,
        absolute_volume_tolerance_m3: 1e-12,
        absolute_amount_tolerance_mol: 1e-12,
        maximum_steps: 50_000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn snapshot(system: &OsmoticNetwork) -> Vec<u64> {
    let mut bits = vec![system.time_s().to_bits()];
    for compartment in system.compartments() {
        bits.push(compartment.solution().solvent_volume_m3().to_bits());
        bits.extend(
            compartment
                .solution()
                .solute_moles()
                .iter()
                .map(|n| n.to_bits()),
        );
    }
    bits
}

fn interacting_pair() -> OsmoticNetwork {
    network(vec![
        cell(
            AqueousSolution::nonideal(model(vec![10.0], vec![2.0]), 1.0, vec![1.0]).unwrap(),
            0.0,
            0.0,
        ),
        cell(ideal(1.0, &[1.0]), 0.0, 0.0),
    ])
}

#[test]
fn g1_zero_interaction_model_recovers_ideal_finite_water_transport() {
    let states = [(1.0, [0.8, 0.2], 0.1 * RT), (2.0, [1.0, 2.0], 0.0)];
    let shared = model(vec![10.0; 2], vec![0.0; 4]);
    let mut reference = network(
        states
            .iter()
            .map(|&(v, n, p)| cell(ideal(v, &n), p, 0.2 * RT))
            .collect(),
    );
    let mut nonideal = network(
        states
            .iter()
            .map(|&(v, n, p)| {
                cell(
                    AqueousSolution::nonideal(shared.clone(), v, n.to_vec()).unwrap(),
                    p,
                    0.2 * RT,
                )
            })
            .collect(),
    );
    reference.advance(0.5, controls(0.05), || false).unwrap();
    let report = nonideal.advance(0.5, controls(0.05), || false).unwrap();
    for (a, b) in reference.compartments().iter().zip(nonideal.compartments()) {
        near(
            a.solution().solvent_volume_m3(),
            b.solution().solvent_volume_m3(),
            2e-10,
        );
        assert_eq!(a.solution().solute_moles(), b.solution().solute_moles());
        near(
            a.solution().water_activity().unwrap(),
            b.solution().water_activity().unwrap(),
            1e-13,
        );
    }
    near(report.water_balance_residual_m3, 0.0, 5e-13);
    assert_eq!(&report.solute_balance_residual_mol[..2], &[0.0, 0.0]);
}

#[test]
fn g1_nonideal_water_reaches_independent_quadratic_equilibrium() {
    let mut system = interacting_pair();
    // nL=nR=1, B_L=2, V_L+V_R=2, no hydrostatic difference:
    // PiL/RT=1/VL+1/VL², PiR/RT=1/(2-VL).
    // Equality reduces to 2 VL² - VL - 2 = 0.
    let equilibrium = (1.0 + 17.0_f64.sqrt()) / 4.0;
    let energy = |system: &OsmoticNetwork| {
        let left = system.compartments()[0].solution().solvent_volume_m3();
        let right = system.compartments()[1].solution().solvent_volume_m3();
        -left.ln() - right.ln() - 2.0 + 1.0 / left
    };
    let mut previous = energy(&system);
    for _ in 0..4 {
        system.advance(0.25, controls(0.1), || false).unwrap();
        let current = energy(&system);
        assert!(current <= previous + 1e-12);
        previous = current;
    }
    assert!(system.compartments()[0].solution().solvent_volume_m3() > 1.0);
    let report = system.advance(99.0, controls(50.0), || false).unwrap();
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    near(left.solvent_volume_m3(), equilibrium, 2e-8);
    near(right.solvent_volume_m3(), 2.0 - equilibrium, 2e-8);
    near(
        left.osmotic_pressure_pa().unwrap(),
        right.osmotic_pressure_pa().unwrap(),
        2e-4,
    );
    near(
        left.water_activity().unwrap(),
        right.water_activity().unwrap(),
        2e-12,
    );
    near(
        left.solvent_volume_m3() + right.solvent_volume_m3(),
        2.0,
        3e-12,
    );
    assert_eq!(left.solute_moles(), &[1.0]);
    assert_eq!(right.solute_moles(), &[1.0]);
    near(report.water_balance_residual_m3, 0.0, 3e-12);
    assert_eq!(report.solute_balance_residual_mol[0], 0.0);
    assert!(energy(&system) <= previous + 1e-12);
}

#[test]
fn g1_cross_interaction_drives_water_despite_equal_ideal_particle_totals() {
    let shared = model(vec![4.0; 2], vec![0.0, 0.2, 0.2, 0.0]);
    let amounts = [[1.0, 1.0], [2.0, 0.0]];
    let mut system = network(
        amounts
            .iter()
            .map(|n| {
                cell(
                    AqueousSolution::nonideal(shared.clone(), 1.0, n.to_vec()).unwrap(),
                    0.0,
                    0.0,
                )
            })
            .collect(),
    );
    let mut ideal_system = network(
        amounts
            .iter()
            .map(|n| cell(ideal(1.0, n), 0.0, 0.0))
            .collect(),
    );
    let dt = 1e-6;
    let report = system.advance(dt, controls(dt), || false).unwrap();
    ideal_system.advance(dt, controls(dt), || false).unwrap();
    // Identical ideal osmotic terms cancel. Only the left c1*c2 cross term
    // remains: (PiL-PiR)/RT=0.2, so dVL/dt=0.2 initially at Lp*RT=1.
    near(
        (system.compartments()[0].solution().solvent_volume_m3() - 1.0) / dt,
        0.2,
        3e-6,
    );
    for (i, initial) in amounts.iter().enumerate() {
        assert_eq!(
            ideal_system.compartments()[i]
                .solution()
                .solvent_volume_m3(),
            1.0
        );
        assert_eq!(system.compartments()[i].solution().solute_moles(), initial);
    }
    near(report.water_balance_residual_m3, 0.0, 3e-14);
    assert_eq!(&report.solute_balance_residual_mol[..2], &[0.0, 0.0]);
}

#[test]
fn g1_very_long_steps_preserve_positive_closed_nonideal_inventories() {
    let mut system = interacting_pair();
    let report = system.advance(1e6, controls(1e6), || false).unwrap();
    let equilibrium = (1.0 + 17.0_f64.sqrt()) / 4.0;
    near(
        system.compartments()[0].solution().solvent_volume_m3(),
        equilibrium,
        2e-7,
    );
    near(
        system
            .compartments()
            .iter()
            .map(|c| c.solution().solvent_volume_m3())
            .sum(),
        2.0,
        3e-12,
    );
    for compartment in system.compartments() {
        assert!(compartment.solution().solvent_volume_m3() > compartment.minimum_volume_m3);
        assert_eq!(compartment.solution().solute_moles(), &[1.0]);
    }
    assert_eq!(system.time_s(), 1e6);
    assert!(report.maximum_scaled_error <= 1.0);
    near(report.water_balance_residual_m3, 0.0, 3e-12);
}

#[test]
fn g0_nonideal_transport_requires_matching_retention_and_film_laws() {
    let nonideal = cell(
        AqueousSolution::nonideal(model(vec![10.0], vec![0.2]), 1.0, vec![1.0]).unwrap(),
        0.0,
        0.0,
    );
    let ideal = cell(ideal(1.0, &[1.0]), 0.0, 0.0);
    for cells in [vec![nonideal.clone(), ideal.clone()], vec![ideal, nonideal]] {
        for pair in [(1.0, 0.1), (0.9, 0.0), (0.0, 0.0)] {
            let link = OsmoticMembrane::new(0, 1, 1.0, law(&[pair])).unwrap();
            assert!(matches!(
                OsmoticNetwork::new(names(1), cells.clone(), vec![link], 1 << 20),
                Err(OsmoticNetworkError::Invalid(_))
            ));
        }
        let film = Some(FilmPermeance::new(0.5).unwrap());
        let material =
            PolarizedMembraneTransport::new(law(&[(1.0, 0.0)]), vec![film], vec![film]).unwrap();
        let link = OsmoticMembrane::new_polarized(0, 1, 1.0, material).unwrap();
        assert!(matches!(
            OsmoticNetwork::new(names(1), cells, vec![link], 1 << 20),
            Err(OsmoticNetworkError::Invalid(_))
        ));
    }
}

#[test]
fn g0_concentration_support_exhaustion_preserves_accepted_state_and_clock() {
    // Left starts exactly at its supported concentration ceiling. The imposed
    // pressure would expel water and immediately exceed that original support.
    let mut system = network(vec![
        cell(
            AqueousSolution::nonideal(model(vec![1.0], vec![2.0]), 1.0, vec![1.0]).unwrap(),
            5.0 * RT,
            0.0,
        ),
        cell(ideal(1.0, &[0.0]), 0.0, 0.0),
    ]);
    let before = snapshot(&system);
    let mut limits = controls(1.0);
    limits.maximum_rejections = 2;
    assert!(system.advance(1.0, limits, || false).is_err());
    assert_eq!(snapshot(&system), before);
    assert!(
        system.compartments()[0]
            .solution()
            .osmotic_pressure_at(0.9, &[1.0])
            .is_err()
    );
}

#[test]
fn g1_unsupported_distant_trial_does_not_hide_a_supported_water_root() {
    let left = cell(
        AqueousSolution::nonideal(model(vec![1.1], vec![2.0]), 1.0, vec![1.0]).unwrap(),
        1.3 * RT,
        0.0,
    );
    let right = cell(ideal(1.0, &[1.0]), 0.0, 0.0);
    // The initial Euler water predictor is 3 m³ over ten seconds. The
    // bracket's capped first trial, ~0.125 m³, exceeds c_max=1.1; the true
    // supported root is smaller than the 1-1/1.1 volume-domain boundary.
    assert!(left.solution().osmotic_pressure_at(0.875, &[1.0]).is_err());
    let endpoint = exchange(
        &left,
        &right,
        &retained(1),
        10.0,
        controls(10.0),
        &mut 0,
        &mut || false,
    )
    .unwrap();
    let water = 1.0 - endpoint.volumes[0];
    assert!(water > 0.0625 && water < 1.0 - 1.0 / 1.1);
    // Eliminating the independently specified rational pressures from the
    // backward-Euler water balance gives w⁴-14w³-8w²+44w-3=0.
    let polynomial = (((water - 14.0) * water - 8.0) * water + 44.0) * water - 3.0;
    near(polynomial, 0.0, 2e-9);
    near(endpoint.volumes[0] + endpoint.volumes[1], 2.0, 1e-14);
    assert_eq!([endpoint.left[0], endpoint.right[0]], [1.0, 1.0]);
}

#[test]
fn g0_late_cancellation_preserves_nonideal_checkpoint() {
    let mut system = interacting_pair();
    system.advance(0.1, controls(0.05), || false).unwrap();
    let before = snapshot(&system);
    let mut polls = 0;
    let result = system.advance(10.0, controls(0.1), || {
        polls += 1;
        polls > 40
    });
    assert!(matches!(result, Err(OsmoticNetworkError::Cancelled)));
    assert!(polls > 40);
    assert_eq!(snapshot(&system), before);
}

#[test]
fn g0_shared_nonideal_model_is_included_in_workspace_admission() {
    let n = 16;
    let mut interaction = vec![0.0; n * n];
    for i in 0..n {
        interaction[i * n + i] = 0.01;
    }
    let shared = model(vec![10.0; n], interaction);
    let amounts = vec![0.1; n];
    let ideal_cells = vec![cell(ideal(1.0, &amounts), 0.0, 0.0); 2];
    let baseline = OsmoticNetwork::new(names(n), ideal_cells, vec![], 1 << 20).unwrap();
    let nonideal_cells = vec![
        cell(
            AqueousSolution::nonideal(shared.clone(), 1.0, amounts).unwrap(),
            0.0,
            0.0
        );
        2
    ];
    let required = match OsmoticNetwork::new(
        names(n),
        nonideal_cells.clone(),
        vec![],
        baseline.workspace_bytes(),
    ) {
        Err(OsmoticNetworkError::Workspace { required_bytes, .. }) => required_bytes,
        other => panic!("shared model storage was omitted: {other:?}"),
    };
    assert!(required > baseline.workspace_bytes() + shared.heap_bytes());
    assert!(matches!(
        OsmoticNetwork::new(names(n), nonideal_cells.clone(), vec![], required - 1),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let accepted = OsmoticNetwork::new(names(n), nonideal_cells, vec![], required).unwrap();
    assert_eq!(accepted.workspace_bytes(), required);
}
