use super::*;
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::{Mesh1d, MixtureGas};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod suspension;

const AREA: f64 = 0.7;
const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-1.0, 1.0],
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 2e-11,
        maximum_iterations: 64,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn left_gas(reference: f64) -> MixtureGas {
    gas(vec![
        source("left", R, 3.5, 0.0, reference),
        source("left-absent", 0.5 * R, 3.0, 0.0, -2.0 * reference),
    ])
}

fn right_gas(reference: f64) -> MixtureGas {
    gas(vec![source("right", 2.0 * R, 4.5, 0.0, 3.0 * reference)])
}

fn chamber(
    gas: MixtureGas,
    start: f64,
    length: f64,
    count: usize,
    area: f64,
    temperature: f64,
) -> MixtureEuler1d {
    let positions: Vec<_> = (0..=count)
        .map(|i| start + length * (i as f64 / count as f64))
        .collect();
    let volumes = positions.windows(2).map(|p| area * (p[1] - p[0])).collect();
    let mesh = Mesh1d::new(positions, vec![area; count + 1], volumes).unwrap();
    let mut densities = vec![0.0; gas.species_count()];
    densities[0] = 1.0;
    let state = gas.conserved(&densities, 0.0, temperature).unwrap();
    MixtureEuler1d::new_on_mesh(gas, mesh, vec![state; count], CLOSED, || false).unwrap()
}

fn pair(
    left_pressure: f64,
    right_pressure: f64,
    reference: f64,
) -> (MixtureEuler1d, MixtureEuler1d) {
    (
        chamber(left_gas(reference), 0.0, 1.0, 12, AREA, left_pressure),
        chamber(right_gas(reference), 1.0, 1.0, 20, AREA, 2.0 * right_pressure),
    )
}

fn step(
    left: &mut MixtureEuler1d,
    right: &mut MixtureEuler1d,
    piston: &mut Piston,
    cancelled: impl FnMut() -> bool,
) -> Result<PistonStep<SharedPistonGasReport>, PistonError> {
    left.step_with_shared_piston(
        right,
        piston,
        Time::new(0.001),
        StepControl::default(),
        0.0,
        controls(),
        cancelled,
    )
}

fn pressure(flow: &MixtureEuler1d, index: usize) -> f64 {
    flow.gas()
        .primitive(flow.cell(index).unwrap(), || false)
        .unwrap()
        .pressure
}

#[test]
fn g1_equal_pressure_different_gases_remain_at_rest_without_species_exchange() {
    let (mut left, mut right) = pair(1.0, 1.0, 0.0);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let before = [left.totals(), right.totals()];
    for _ in 0..4 {
        let report = step(&mut left, &mut right, &mut piston, || false).unwrap();
        assert!(report.admitted(controls()));
        close(piston.velocity_m_s(), 0.0, 2e-13);
        close(report.displacement_m, 0.0, 2e-15);
        close(report.gas_impulse_n_s, 0.0, 2e-15);
        close(report.fixed_wall_impulse_n_s, 0.0, 2e-15);
        assert_eq!(report.exterior_pressure_work_j, 0.0);
        let dt = report.gas.left.elapsed.value();
        // Independent p*A*dt on opposite faces of the impermeable piston.
        close(report.gas.left.right_boundary_change.momentum, -AREA * dt, 2e-15);
        close(report.gas.right.left_boundary_change.momentum, AREA * dt, 2e-15);
        assert_eq!(report.gas.left.elapsed, report.gas.right.elapsed);
        assert_eq!(left.time(), right.time());
        for (side, original) in [(&left, before[0]), (&right, before[1])] {
            let after = side.totals();
            for k in 0..side.gas().species_count() {
                close(after.species_masses[k], original.species_masses[k], 2e-14);
            }
            close(after.energy, original.energy, 5e-14);
            for i in 0..side.cell_count() {
                close(pressure(side, i), 1.0, 2e-13);
            }
        }
        for i in 0..left.cell_count() {
            assert_eq!(left.cell(i).unwrap().partial_densities[1], 0.0);
        }
    }
}

#[test]
fn g1_both_resolved_pressures_feed_back_into_motion_work_and_support_balance() {
    let (mut left, mut right) = pair(2.0, 1.0, 0.0);
    let (mut shifted_left, mut shifted_right) = pair(2.0, 1.0, 1e18);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut shifted_piston = piston;
    let before = [left.totals(), right.totals()];
    let mut supports = 0.0;
    let mut work = 0.0;
    let mut imposed_impulse = 0.0;
    let external = -0.03;
    for _ in 0..40 {
        let report = left
            .step_with_shared_piston(
                &mut right, &mut piston, Time::new(0.001), StepControl::default(),
                external, controls(), || false,
            )
            .unwrap();
        let shifted = shifted_left
            .step_with_shared_piston(
                &mut shifted_right, &mut shifted_piston, Time::new(0.001),
                StepControl::default(), external, controls(), || false,
            )
            .unwrap();
        assert_eq!(report, shifted, "separate source energy datums must not drive the piston");
        assert_eq!(left.states, shifted_left.states);
        assert_eq!(right.states, shifted_right.states);
        assert_eq!(piston, shifted_piston);
        assert!(report.admitted(controls()));
        let gas = report.gas;
        assert!(report.midpoint_velocity_m_s > 0.0);
        assert!(gas.left.right_boundary_change.energy < 0.0);
        assert!(gas.right.left_boundary_change.energy > 0.0);
        assert_eq!(gas.left.left_boundary_change.energy, 0.0);
        assert_eq!(gas.right.right_boundary_change.energy, 0.0);
        for side in [gas.left, gas.right] {
            close(side.conservation_residual.energy, 0.0, 2e-13);
            for mass in side.boundary_change.species_masses {
                assert_eq!(mass, 0.0);
            }
        }
        close(
            gas.left.right_boundary_change.energy + gas.right.left_boundary_change.energy,
            -report.midpoint_velocity_m_s * report.gas_impulse_n_s,
            2e-14,
        );
        assert_eq!(left.mesh().positions()[left.cell_count()], right.mesh().positions()[0]);
        close(
            left.mesh().volumes().iter().sum::<f64>()
                + right.mesh().volumes().iter().sum::<f64>(),
            2.0 * AREA,
            3e-14,
        );
        supports += report.fixed_wall_impulse_n_s;
        work += report.external_force_work_j;
        imposed_impulse += external * gas.left.elapsed.value();
    }
    let after = [left.totals(), right.totals()];
    for (side, initial, final_state) in [
        (&left, before[0], after[0]),
        (&right, before[1], after[1]),
    ] {
        for k in 0..side.gas().species_count() {
            close(final_state.species_masses[k], initial.species_masses[k], 3e-13);
        }
    }
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    close(
        after[0].energy - before[0].energy + after[1].energy - before[1].energy
            + kinetic - work,
        0.0,
        8e-10,
    );
    close(
        after[0].momentum - before[0].momentum + after[1].momentum - before[1].momentum
            + piston.mass_kg() * piston.velocity_m_s() + supports - imposed_impulse,
        0.0,
        8e-11,
    );
    assert!(piston.velocity_m_s() > 0.02);
    assert!(pressure(&left, left.cell_count() - 1) < pressure(&left, 0));
    assert!(pressure(&right, 0) > pressure(&right, right.cell_count() - 1));
    assert_eq!(left.mesh().positions()[0], 0.0);
    assert_eq!(right.mesh().positions()[right.cell_count()], 2.0);

    // Reversing the pressure difference drives the same two-sided wall left.
    let (mut left, mut right) = pair(1.0, 2.0, 0.0);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    for _ in 0..4 {
        let report = step(&mut left, &mut right, &mut piston, || false).unwrap();
        assert!(report.displacement_m < 0.0);
        assert!(report.gas.left.right_boundary_change.energy > 0.0);
        assert!(report.gas.right.left_boundary_change.energy < 0.0);
        assert!(report.admitted(controls()));
    }
    assert!(piston.velocity_m_s() < 0.0);
}

#[test]
fn g1_common_right_chamber_cfl_shortening_replays_both_original_moving_fields() {
    let make = || (
        chamber(left_gas(0.0), 0.0, 1.0, 8, AREA, 2.0),
        chamber(right_gas(0.0), 1.0, 0.1, 40, AREA, 1.6),
    );
    let (mut left, mut right) = make();
    let (mut replay_left, mut replay_right) = make();
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let report = left
        .step_with_shared_piston(
            &mut right, &mut piston, Time::new(0.1), StepControl::default(), 0.0,
            controls(), || false,
        )
        .unwrap();
    assert!(report.gas.left.elapsed.value() < 0.001);
    assert_eq!(report.gas.left.elapsed, report.gas.right.elapsed);
    assert_eq!(report.gas.left.halvings, report.gas.right.halvings);
    assert!(report.gas.left.maximum_cfl <= 0.40000000000001);
    assert!(report.gas.right.maximum_cfl <= 0.40000000000001);
    for (reference, accepted, side, speeds) in [
        (&mut replay_left, &left, report.gas.left, [0.0, report.midpoint_velocity_m_s]),
        (&mut replay_right, &right, report.gas.right, [report.midpoint_velocity_m_s, 0.0]),
    ] {
        let replay = reference
            .step_with_moving_walls(
                side.elapsed,
                StepControl {
                    reconstruction: if side.first_order_retry {
                        Reconstruction::FirstOrder
                    } else {
                        Reconstruction::Muscl
                    },
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                speeds,
                || false,
            )
            .unwrap();
        assert_eq!(replay.elapsed, side.elapsed);
        assert_eq!(reference.mesh(), accepted.mesh());
        assert_eq!(reference.time(), accepted.time());
        for (a, b) in reference.states.iter().zip(&accepted.states) {
            close(*a, *b, 4e-14);
        }
        close(replay.left_boundary_change.energy, side.left_boundary_change.energy, 2e-14);
        close(replay.right_boundary_change.energy, side.right_boundary_change.energy, 2e-14);
    }
    assert_eq!(left.time(), right.time());
}

#[test]
fn g4_late_cancel_and_coupling_refusals_preserve_both_chambers_and_replay() {
    let (mut left, mut right) = pair(2.0, 1.0, 0.0);
    let (mut reference_left, mut reference_right) = pair(2.0, 1.0, 0.0);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut reference_piston = piston;
    step(&mut left, &mut right, &mut piston, || false).unwrap();
    step(&mut reference_left, &mut reference_right, &mut reference_piston, || false).unwrap();
    let left_state = left.states.clone();
    let right_state = right.states.clone();
    let left_mesh = left.mesh().clone();
    let right_mesh = right.mesh().clone();
    let time = left.time();
    let original_piston = piston;
    let mut polls = 0;
    let expected = step(&mut reference_left, &mut reference_right, &mut reference_piston, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 100);
    let mut count = 0;
    let failure = step(&mut left, &mut right, &mut piston, || {
        count += 1;
        count == polls
    });
    assert_eq!(failure, Err(PistonError::Cancelled));
    assert_eq!(count, polls);
    assert_eq!(left.states, left_state);
    assert_eq!(right.states, right_state);
    assert_eq!(left.mesh(), &left_mesh);
    assert_eq!(right.mesh(), &right_mesh);
    assert_eq!(left.time(), time);
    assert_eq!(right.time(), time);
    assert_eq!(piston, original_piston);
    let actual = step(&mut left, &mut right, &mut piston, || false).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(left.states, reference_left.states);
    assert_eq!(right.states, reference_right.states);
    assert_eq!(left.mesh(), reference_left.mesh());
    assert_eq!(right.mesh(), reference_right.mesh());
    assert_eq!(piston, reference_piston);

    for refusal in 0..4 {
        let (mut left, mut right) = pair(2.0, 1.0, 0.0);
        match refusal {
            0 => right = chamber(right_gas(0.0), 1.01, 1.0, 20, AREA, 2.0),
            1 => right = chamber(right_gas(0.0), 1.0, 1.0, 20, 2.0 * AREA, 2.0),
            2 => right.time = 0.01,
            _ => right.boundaries = Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Transmissive,
            },
        }
        let left_state = left.states.clone();
        let right_state = right.states.clone();
        let times = [left.time(), right.time()];
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        assert!(matches!(
            step(&mut left, &mut right, &mut piston, || false),
            Err(PistonError::Invalid(_))
        ));
        assert_eq!(left.states, left_state);
        assert_eq!(right.states, right_state);
        assert_eq!([left.time(), right.time()], times);
        assert_eq!(piston.velocity_m_s(), 0.0);
    }
}
