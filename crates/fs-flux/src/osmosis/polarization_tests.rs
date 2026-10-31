//! Independent film-limited transport, finite-inventory, and refusal regressions.

use super::*;
use fs_material::osmosis::polarization::FilmPermeance;
use fs_material::osmosis::{HydraulicPermeance, ReflectionCoefficient, SolutePermeance};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const T: f64 = 298.15;
const RT: f64 = R * T;

fn cell(volume: f64, amounts: &[f64], pressure: f64, elastance: f64) -> OsmoticCompartment {
    OsmoticCompartment::new(
        IdealDiluteSolution::new(T, 18e-6, volume, amounts.to_vec()).unwrap(),
        OsmoticPressureVolume::new(volume, pressure, elastance).unwrap(),
        volume * 1e-6,
    )
    .unwrap()
}

fn material(lp: f64, species: &[(f64, f64)]) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        species
            .iter()
            .map(|&(s, p)| {
                (
                    ReflectionCoefficient::new(s).unwrap(),
                    SolutePermeance::new(p).unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
}

fn polarized(
    left: usize,
    right: usize,
    law: MembraneTransport,
    kl: &[Option<f64>],
    kr: &[Option<f64>],
) -> OsmoticMembrane {
    let prepare = |values: &[Option<f64>]| {
        values
            .iter()
            .map(|k| k.map(|k| FilmPermeance::new(k).unwrap()))
            .collect()
    };
    OsmoticMembrane::new_polarized(
        left,
        right,
        1.0,
        PolarizedMembraneTransport::new(law, prepare(kl), prepare(kr)).unwrap(),
    )
    .unwrap()
}

fn network(cells: Vec<OsmoticCompartment>, membranes: Vec<OsmoticMembrane>) -> OsmoticNetwork {
    let names = (0..cells[0].solution().solute_moles().len())
        .map(|i| format!("solute-{i}"))
        .collect();
    OsmoticNetwork::new(names, cells, membranes, 1 << 20).unwrap()
}

fn controls(step: f64) -> OsmoticStepControls {
    OsmoticStepControls {
        maximum_step_s: step,
        relative_tolerance: 1e-6,
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
    for cell in system.compartments() {
        bits.push(cell.solution().solvent_volume_m3().to_bits());
        bits.extend(cell.solution().solute_moles().iter().map(|x| x.to_bits()));
    }
    bits
}

#[test]
fn absent_films_preserve_the_original_network_bit_for_bit() {
    let cells = vec![
        cell(1.0, &[1.0, 0.0], 0.4 * RT, 0.5 * RT),
        cell(1.5, &[3.0, 2.0], 0.0, 0.3 * RT),
    ];
    let law = material(0.2 / RT, &[(0.8, 0.1), (0.2, 0.4)]);
    let mut bare = network(
        cells.clone(),
        vec![OsmoticMembrane::new(0, 1, 1.0, law.clone()).unwrap()],
    );
    let mut prepared = network(
        cells,
        vec![polarized(0, 1, law, &[None, None], &[None, None])],
    );
    let a = bare.advance(0.3, controls(0.02), || false).unwrap();
    let b = prepared.advance(0.3, controls(0.02), || false).unwrap();
    assert_eq!(snapshot(&bare), snapshot(&prepared));
    assert_eq!(
        (a.accepted_steps, a.rejected_steps),
        (b.accepted_steps, b.rejected_steps)
    );
    assert_eq!(a.residual_evaluations, b.residual_evaluations);
}

#[test]
fn zero_water_flow_matches_series_diffusion_into_an_empty_receiver() {
    let mut system = network(
        vec![cell(1.0, &[30.0], 0.0, 0.0), cell(2.0, &[0.0], 0.0, 0.0)],
        vec![polarized(
            0,
            1,
            material(0.0, &[(0.7, 0.5)]),
            &[Some(0.25)],
            &[Some(1.0)],
        )],
    );
    let report = system.advance(4.0, controls(0.1), || false).unwrap();
    let effective = 1.0 / (1.0 / 0.25 + 1.0 / 0.5 + 1.0);
    let difference = 30.0 * (-effective * 1.5_f64 * 4.0).exp();
    let left = system.compartments()[0].solution();
    let right = system.compartments()[1].solution();
    near(
        left.solute_moles()[0],
        10.0 + (2.0 / 3.0) * difference,
        0.01,
    );
    near(
        left.solute_moles()[0] + right.solute_moles()[0],
        30.0,
        3e-12,
    );
    assert!(right.solute_moles()[0] > 0.0);
    assert_eq!(
        (left.solvent_volume_m3(), right.solvent_volume_m3()),
        (1.0, 2.0)
    );
    assert_eq!(report.water_balance_residual_m3, 0.0);
    near(report.solute_balance_residual_mol[0], 0.0, 3e-12);
}

#[test]
fn tiny_finite_step_matches_independent_polarization_limited_water_rates() {
    let dt = 1e-7;
    for (lp_rt, expected_rate) in [
        (1.0, 0.723_654_785_457_764_8),
        (10.0, 0.900_237_151_872_033_6),
        (1e6, 0.924_195_996_234_962_6),
    ] {
        let mut system = network(
            vec![cell(1.0, &[1.0], 0.0, 0.0), cell(1.0, &[4.0], 0.0, 0.0)],
            vec![polarized(
                0,
                1,
                material(lp_rt / RT, &[(1.0, 0.0)]),
                &[Some(1.0)],
                &[Some(2.0)],
            )],
        );
        let report = system.advance(dt, controls(dt), || false).unwrap();
        let left = system.compartments()[0].solution();
        let right = system.compartments()[1].solution();
        // References solve J=q*(4 exp(-J/2)-exp(J)) for fixed reservoirs.
        // Finite reservoirs shift the interval rate by O(dt), not exactly zero.
        near((right.solvent_volume_m3() - 1.0) / dt, expected_rate, 2e-6);
        near(
            left.solvent_volume_m3() + right.solvent_volume_m3(),
            2.0,
            1e-14,
        );
        assert_eq!(left.solute_moles(), &[1.0]);
        assert_eq!(right.solute_moles(), &[4.0]);
        assert_eq!(report.accepted_steps, 1);
    }
}

#[test]
fn films_slow_the_transient_without_changing_retained_solute_equilibrium() {
    let cells = vec![cell(1.0, &[1.0], 0.0, 0.0), cell(1.0, &[2.0], 0.0, 0.0)];
    let law = material(0.5 / RT, &[(1.0, 0.0)]);
    let mut bare = network(
        cells.clone(),
        vec![OsmoticMembrane::new(0, 1, 1.0, law.clone()).unwrap()],
    );
    let mut limited = network(
        cells,
        vec![polarized(0, 1, law, &[Some(0.1)], &[Some(0.2)])],
    );
    bare.advance(0.5, controls(0.1), || false).unwrap();
    limited.advance(0.5, controls(0.1), || false).unwrap();
    let bare_transfer = bare.compartments()[1].solution().solvent_volume_m3() - 1.0;
    let film_transfer = limited.compartments()[1].solution().solvent_volume_m3() - 1.0;
    assert!(film_transfer > 0.0 && film_transfer < 0.5 * bare_transfer);
    for system in [&mut bare, &mut limited] {
        system.advance(120.0, controls(1.0), || false).unwrap();
        let left = system.compartments()[0].solution();
        let right = system.compartments()[1].solution();
        near(left.solvent_volume_m3(), 2.0 / 3.0, 2e-7);
        near(right.solvent_volume_m3(), 4.0 / 3.0, 2e-7);
        near(
            left.solvent_volume_m3() + right.solvent_volume_m3(),
            2.0,
            3e-12,
        );
        assert_eq!(left.solute_moles(), &[1.0]);
        assert_eq!(right.solute_moles(), &[2.0]);
    }
}

fn free_energy(system: &OsmoticNetwork) -> f64 {
    system
        .compartments()
        .iter()
        .map(|cell| {
            let v = cell.solution().solvent_volume_m3();
            let dv = v - cell.pressure.reference_volume_m3;
            let mechanical = cell.pressure.reference_pressure_pa * dv
                + 0.5 * cell.pressure.elastance_pa_per_m3 * dv * dv;
            mechanical / RT
                + cell
                    .solution()
                    .solute_moles()
                    .iter()
                    .map(|&n| {
                        if n == 0.0 {
                            0.0
                        } else {
                            n * ((n / v).ln() - 1.0)
                        }
                    })
                    .sum::<f64>()
        })
        .sum()
}

#[test]
fn leaking_multispecies_films_conserve_inventories_dissipate_and_reverse_orientation() {
    let cells = vec![
        cell(1.0, &[2.0, 0.0, 1.0], 0.4 * RT, 0.5 * RT),
        cell(1.0, &[0.5, 3.0, 0.2], 0.0, 0.5 * RT),
    ];
    let law = material(0.2 / RT, &[(0.8, 0.1), (0.2, 0.4), (1.0, 0.02)]);
    let kl = [Some(0.3), None, Some(0.8)];
    let kr = [Some(0.6), Some(0.7), None];
    let mut forward = network(cells.clone(), vec![polarized(0, 1, law.clone(), &kl, &kr)]);
    let mut reverse = network(cells, vec![polarized(1, 0, law, &kr, &kl)]);
    let mut previous = free_energy(&forward);
    for _ in 0..8 {
        let report = forward.advance(0.25, controls(0.05), || false).unwrap();
        reverse.advance(0.25, controls(0.05), || false).unwrap();
        let energy = free_energy(&forward);
        assert!(energy <= previous + 1e-10 * (1.0 + previous.abs()));
        previous = energy;
        near(
            forward
                .compartments()
                .iter()
                .map(|c| c.solution().solvent_volume_m3())
                .sum(),
            2.0,
            3e-12,
        );
        for (i, total) in [2.5, 3.0, 1.2].into_iter().enumerate() {
            near(
                forward
                    .compartments()
                    .iter()
                    .map(|c| c.solution().solute_moles()[i])
                    .sum(),
                total,
                3e-12,
            );
            near(report.solute_balance_residual_mol[i], 0.0, 3e-12);
        }
        for (a, b) in forward.compartments().iter().zip(reverse.compartments()) {
            near(
                a.solution().solvent_volume_m3(),
                b.solution().solvent_volume_m3(),
                5e-10,
            );
            for (&na, &nb) in a
                .solution()
                .solute_moles()
                .iter()
                .zip(b.solution().solute_moles())
            {
                assert!(na >= 0.0 && nb >= 0.0);
                near(na, nb, 5e-10);
            }
        }
    }
    assert!(forward.compartments()[0].solution().solute_moles()[1] > 0.1);
}

#[test]
fn polarization_heap_storage_and_late_cancellation_respect_atomic_bounds() {
    let cells = vec![cell(1.0, &[3.0], 0.0, 0.0), cell(1.0, &[1.0], 0.0, 0.0)];
    let law = material(0.0, &[(0.5, 0.2)]);
    let bare = network(
        cells.clone(),
        vec![OsmoticMembrane::new(0, 1, 1.0, law.clone()).unwrap()],
    );
    let links = vec![polarized(0, 1, law, &[Some(0.3)], &[Some(0.5)])];
    let mut system = network(cells.clone(), links.clone());
    assert!(
        system.workspace_bytes()
            >= bare.workspace_bytes() + 2 * core::mem::size_of::<Option<FilmPermeance>>()
    );
    let refused = OsmoticNetwork::new(
        vec!["solute-0".into()],
        cells,
        links,
        bare.workspace_bytes(),
    );
    assert!(matches!(
        refused,
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let mut control = controls(0.1);
    control.relative_tolerance = 0.9;
    let before = snapshot(&system);
    let mut uninterrupted = system.clone();
    let mut polls = 0;
    let result = system.advance(1.0, control, || {
        polls += 1;
        polls == 10
    });
    assert!(matches!(result, Err(OsmoticNetworkError::Cancelled)));
    assert_eq!(snapshot(&system), before);
    // Jv=0 needs one residual per exchange; poll10 occurs after an accepted
    // first temporal step, so rollback covers already staged physical progress.
    assert_ne!(
        system.work[0].solution().solute_moles(),
        system.compartments()[0].solution().solute_moles()
    );
    system.advance(1.0, control, || false).unwrap();
    uninterrupted.advance(1.0, control, || false).unwrap();
    assert_eq!(snapshot(&system), snapshot(&uninterrupted));
}

#[test]
fn linear_compliance_is_resolved_with_a_small_nonlinear_iteration_budget() {
    let left = cell(1.0, &[], 1.0, 2.0);
    let right = cell(1.0, &[], 0.0, 3.0);
    let link = polarized(0, 1, material(1.0, &[]), &[], &[]);
    let mut control = controls(0.2);
    control.maximum_root_iterations = 4;
    let mut evaluations = 0;
    let endpoint = exchange(
        &left,
        &right,
        &link,
        0.2,
        control,
        &mut evaluations,
        &mut || false,
    )
    .unwrap();
    // Exact backward Euler x=dt*Lp*dp/(1+dt*Lp*(kL+kR)).
    near(endpoint.volumes[0], 0.9, 2e-14);
    near(endpoint.volumes[1], 1.1, 2e-14);
    assert!(evaluations <= 6);
}
