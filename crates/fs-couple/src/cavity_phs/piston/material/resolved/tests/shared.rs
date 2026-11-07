//! Source inertia, two gas fields and one finite thermal store stay coherent.
use super::*;
use fs_gas::multicomponent::{
    Mesh1d, SharedPistonThermal, SharedPistonThermalError,
    wall_exchange::{SharedWallControl, SharedWallStorage, WallControl, WallThermal},
};

fn pair() -> (MixtureEuler1d, MixtureEuler1d) {
    let seed = mixture();
    let chamber = |start: f64, count: usize, cell: usize, pressure: f64| {
        let gas = seed.gas().clone();
        let state = gas.conserved(seed.cell(cell).unwrap().partial_densities, 0.0, pressure).unwrap();
        let positions: Vec<_> = (0..=count)
            .map(|i| start + 0.08 * (i as f64 / count as f64)).collect();
        let volumes = positions.windows(2).map(|x| 0.01 * (x[1] - x[0])).collect();
        let mesh = Mesh1d::new(positions, vec![0.01; count + 1], volumes).unwrap();
        MixtureEuler1d::new_on_mesh(gas, mesh, vec![state; count], CLOSED, || false).unwrap()
    };
    (chamber(0.0, 2, 0, 100000.0), chamber(0.08, 3, 7, 80000.0))
}

fn finite(left: &MixtureEuler1d, right: &MixtureEuler1d) -> SharedPistonThermal<'static> {
    SharedPistonThermal::new(
        left, right,
        SharedWallStorage {
            thermal: WallThermal::Lumped {
                temperature_k: 600.0, heat_capacity_j_per_k: 0.6,
                conductance_w_per_k: 0.0,
            },
            heat_storage: None,
        },
        [20.0, 5.0], SharedWallControl::default(), || false,
    ).unwrap()
}

fn heat_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-12,
        ..WallControl::default()
    }
}

fn temperature(wall: &SharedPistonThermal<'_>) -> f64 {
    match wall.storage().thermal {
        WallThermal::Lumped { temperature_k, .. } => temperature_k,
        _ => unreachable!("the fixture retains its one finite capacity"),
    }
}

fn cells(flow: &MixtureEuler1d) -> Vec<(Vec<f64>, f64, f64)> {
    (0..flow.cell_count()).map(|i| {
        let q = flow.cell(i).unwrap();
        (q.partial_densities.to_vec(), q.momentum, q.energy)
    }).collect()
}

fn same_field(actual: &MixtureEuler1d, expected: &MixtureEuler1d) {
    assert_eq!(cells(actual), cells(expected));
    assert_eq!(actual.mesh(), expected.mesh());
    assert_eq!(actual.time(), expected.time());
}

fn position(flow: &MixtureEuler1d) -> f64 {
    *flow.mesh().positions().last().unwrap()
}

fn energy(
    left: &MixtureEuler1d, right: &MixtureEuler1d, heat: &SharedPistonThermal<'_>,
    wall: &MaterialResolvedPiston, suspended: bool,
) -> f64 {
    left.totals().energy + right.totals().energy + 0.6 * temperature(heat)
        + 0.5 * wall.body().inertia().mass() * wall.piston().velocity_m_s().powi(2)
        + if suspended { 0.5 * 8000.0 * (position(left) - 0.079).powi(2) } else { 0.0 }
}

#[test]
fn both_resolved_pressures_use_one_sourced_inertia_and_original_species() {
    for suspended in [false, true] {
        let specimen = body(4000.0);
        let mut wall = MaterialResolvedPiston::new(specimen.clone(), 0.08).unwrap();
        let mut native = *wall.piston();
        let (mut left, mut right) = pair();
        let (mut a, mut b) = pair();
        let before = [left.totals(), right.totals()];
        for _ in 0..3 {
            let report = if suspended {
                let actual = wall.step_shared_mixture_suspended(
                    &mut left, &mut right, DT, StepControl::default(), 3.0,
                    suspension(), controls(), || false,
                ).unwrap();
                let expected = a.step_with_shared_suspended_piston(
                    &mut b, &mut native, DT, StepControl::default(), 3.0,
                    suspension(), controls(), || false,
                ).unwrap();
                assert_eq!(actual, expected);
                actual.piston
            } else {
                let actual = wall.step_shared_mixture(
                    &mut left, &mut right, DT, StepControl::default(), 3.0,
                    controls(), || false,
                ).unwrap();
                let expected = a.step_with_shared_piston(
                    &mut b, &mut native, DT, StepControl::default(), 3.0,
                    controls(), || false,
                ).unwrap();
                assert_eq!(actual, expected);
                actual
            };
            assert_eq!(wall.body(), &specimen);
            assert_eq!(wall.piston(), &native);
            same_field(&left, &a);
            same_field(&right, &b);
            assert_eq!(position(&left), right.mesh().positions()[0]);
            assert_eq!(left.time(), right.time());
            assert_eq!(report.exterior_pressure_work_j, 0.0);
            assert!(report.gas.left.right_boundary_change.energy < 0.0);
            assert!(report.gas.right.left_boundary_change.energy > 0.0);
            for (flow, initial) in [(&left, &before[0]), (&right, &before[1])] {
                for (after, start) in flow.totals().species_masses.iter().zip(&initial.species_masses) {
                    assert!((after - start).abs() < 1e-14);
                }
            }
        }
    }
}

#[test]
fn finite_wall_heat_and_returned_damping_reach_both_original_pressure_fields() {
    for suspended in [false, true] {
        let specimen = body(4000.0);
        let mut wall = MaterialResolvedPiston::new(specimen.clone(), 0.08).unwrap();
        let mut native = *wall.piston();
        let (mut left, mut right) = pair();
        let (mut a, mut b) = pair();
        let (mut heat, mut reference_heat) = (finite(&left, &right), finite(&a, &b));
        for _ in 0..3 {
            let initial_energy = energy(&left, &right, &heat, &wall, suspended);
            let x0 = position(&left);
            let (report, loss) = if suspended {
                let actual = wall.step_shared_mixture_suspended_with_heat(
                    &mut left, &mut right, &mut heat, DT, StepControl::default(),
                    heat_controls(), 3.0, suspension(), 0.35, controls(), || false,
                ).unwrap();
                let expected = a.step_with_shared_suspended_piston_and_heat(
                    &mut b, &mut reference_heat, &mut native, DT, StepControl::default(),
                    heat_controls(), 3.0, suspension(), 0.35, controls(), || false,
                ).unwrap();
                assert_eq!(actual, expected);
                assert!(actual.damping_energy_j > 0.0);
                assert!((actual.piston.gas.heat.damping_heat_to_wall_j
                    - 0.35 * actual.damping_energy_j).abs() < 1e-14);
                (actual.piston, 0.65 * actual.damping_energy_j)
            } else {
                let actual = wall.step_shared_mixture_with_heat(
                    &mut left, &mut right, &mut heat, DT, StepControl::default(),
                    heat_controls(), 3.0, controls(), || false,
                ).unwrap();
                let expected = a.step_with_shared_piston_and_heat(
                    &mut b, &mut reference_heat, &mut native, DT, StepControl::default(),
                    heat_controls(), 3.0, controls(), || false,
                ).unwrap();
                assert_eq!(actual, expected);
                (actual, 0.0)
            };
            assert_eq!(wall.body(), &specimen);
            assert_eq!(wall.piston(), &native);
            same_field(&left, &a);
            same_field(&right, &b);
            assert_eq!(heat.storage(), reference_heat.storage());
            assert_eq!(heat.last_report(), reference_heat.last_report());
            assert!(report.gas.heat.heat_to_gas_j.iter().all(|q| *q > 0.0));
            let balance = energy(&left, &right, &heat, &wall, suspended)
                - initial_energy + loss - 3.0 * (position(&left) - x0);
            assert!(balance.abs() <= controls().energy_tolerance_j, "{balance:e}");
        }
    }
}

#[test]
fn late_shared_heat_cancellation_preserves_every_owner_and_retries_exactly() {
    let mut wall = MaterialResolvedPiston::new(body(4000.0), 0.08).unwrap();
    let mut native = *wall.piston();
    let (mut left, mut right) = pair();
    let (mut a, mut b) = pair();
    let (mut heat, mut reference_heat) = (finite(&left, &right), finite(&a, &b));
    wall.step_shared_mixture_suspended_with_heat(
        &mut left, &mut right, &mut heat, DT, StepControl::default(),
        heat_controls(), 3.0, suspension(), 0.35, controls(), || false,
    ).unwrap();
    a.step_with_shared_suspended_piston_and_heat(
        &mut b, &mut reference_heat, &mut native, DT, StepControl::default(),
        heat_controls(), 3.0, suspension(), 0.35, controls(), || false,
    ).unwrap();
    let original = wall.clone();
    let states = [cells(&left), cells(&right)];
    let meshes = [left.mesh().clone(), right.mesh().clone()];
    let times = [left.time(), right.time()];
    let storage = heat.storage().clone();
    let last_report = heat.last_report();
    let mut polls = 0;
    let expected = a.step_with_shared_suspended_piston_and_heat(
        &mut b, &mut reference_heat, &mut native, DT, StepControl::default(),
        heat_controls(), 3.0, suspension(), 0.35, controls(),
        || { polls += 1; false },
    ).unwrap();
    let mut visits = 0;
    let error = wall.step_shared_mixture_suspended_with_heat(
        &mut left, &mut right, &mut heat, DT, StepControl::default(),
        heat_controls(), 3.0, suspension(), 0.35, controls(),
        || { visits += 1; visits == polls },
    ).unwrap_err();
    assert_eq!(error, SharedPistonThermalError::Cancelled);
    assert_eq!(wall, original);
    assert_eq!([cells(&left), cells(&right)], states);
    assert_eq!([left.mesh(), right.mesh()], [&meshes[0], &meshes[1]]);
    assert_eq!([left.time(), right.time()], times);
    assert_eq!(heat.storage(), &storage);
    assert_eq!(heat.last_report(), last_report);
    let actual = wall.step_shared_mixture_suspended_with_heat(
        &mut left, &mut right, &mut heat, DT, StepControl::default(),
        heat_controls(), 3.0, suspension(), 0.35, controls(), || false,
    ).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(wall.body(), original.body());
    assert_eq!(wall.piston(), &native);
    same_field(&left, &a);
    same_field(&right, &b);
    assert_eq!(heat.storage(), reference_heat.storage());
    assert_eq!(heat.last_report(), reference_heat.last_report());
}
