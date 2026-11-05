//! Compression, independent chemical gauges and one finite wall inventory.

use super::*;
use crate::dynamics::piston::{PistonStep, PistonSuspension};
use crate::multicomponent::wall_exchange::{
    SharedWallControl, SharedWallStorage, WallControl, WallThermal,
};
use crate::multicomponent::{SharedPistonThermal, SharedPistonThermalStepReport};

const WALL_CAPACITY: f64 = 0.6;

fn heat_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-12,
        relative_energy_tolerance: 1e-12,
        ..WallControl::default()
    }
}

fn wall(
    left: &ReactingEuler1d,
    right: &ReactingEuler1d,
    temperature: f64,
    conductances: [f64; 2],
) -> SharedPistonThermal<'static> {
    SharedPistonThermal::new(
        left.flow(),
        right.flow(),
        SharedWallStorage {
            thermal: WallThermal::Lumped {
                temperature_k: temperature,
                heat_capacity_j_per_k: WALL_CAPACITY,
                conductance_w_per_k: 0.0,
            },
            heat_storage: None,
        },
        conductances,
        SharedWallControl::default(),
        || false,
    )
    .unwrap()
}

fn wall_temperature(wall: &SharedPistonThermal<'_>) -> f64 {
    let WallThermal::Lumped { temperature_k, .. } = wall.storage().thermal else {
        panic!("test wall has original constant capacity");
    };
    temperature_k
}

fn advance(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    wall: &mut SharedPistonThermal<'_>,
    piston: &mut Piston,
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveSharedPistonThermalStep, ReactiveError> {
    left.step_with_shared_piston_and_heat(
        right,
        wall,
        piston,
        Time::new(duration),
        ReactiveStepControl::default(),
        heat_controls(),
        0.01,
        piston_controls(),
        cancelled,
    )
}

fn compare_fields(actual: &ReactingEuler1d, expected: &ReactingEuler1d, tolerance: f64) {
    assert_eq!(actual.flow().time(), expected.flow().time());
    for (&a, &b) in actual.accepted.states.iter().zip(&expected.accepted.states) {
        close(a, b, tolerance);
    }
    for (&a, &b) in actual
        .flow()
        .mesh()
        .volumes()
        .iter()
        .zip(expected.flow().mesh().volumes())
    {
        close(a, b, 2e-11);
    }
    for (&a, &b) in actual
        .flow()
        .mesh()
        .positions()
        .iter()
        .zip(expected.flow().mesh().positions())
    {
        close(a, b, 2e-11);
    }
}

#[test]
fn g1_shared_reacting_thermal_compression_closes_each_source_and_one_wall() {
    let (mut left, mut right) = pair();
    let mut thermal = wall(&left, &right, 4.0, [0.9, 0.3]);
    let mut piston = Piston::new(0.2, 0.04).unwrap();
    let before = [left.flow().totals(), right.flow().totals()];
    let initial_kinetic = 0.5 * 0.2 * 0.04_f64.powi(2);
    let mut drive = 0.0;
    for _ in 0..12 {
        let report = advance(
            &mut left,
            &mut right,
            &mut thermal,
            &mut piston,
            0.002,
            || false,
        )
        .unwrap();
        let middle = report.piston.gas;
        assert_eq!(report.elapsed, middle.transport.left.elapsed);
        assert_eq!(report.elapsed, middle.transport.right.elapsed);
        assert_eq!(left.flow().time(), right.flow().time());
        assert_eq!(thermal.last_report(), Some(middle.heat));
        for stage in middle.heat.stages {
            assert_eq!(stage.elapsed.value(), 0.5 * report.elapsed.value());
            assert_eq!(stage.mechanical_work_to_gas_j, 0.0);
            assert_eq!(stage.external_heat_to_gas_j, 0.0);
        }
        for (j, residual) in [
            report.conservation_residual.left,
            report.conservation_residual.right,
        ]
        .into_iter()
        .enumerate()
        {
            assert_eq!(report.chemistry[j].cell_solves, [16, 24][j]);
            assert!(report.chemistry[j].referenced_energy_change_j > 0.0);
            close(report.mass_residual_kg[j], 0.0, 3e-12);
            for &change in &residual.species_masses[..residual.species_count] {
                close(change, 0.0, 3e-12);
            }
            close(residual.energy, 0.0, piston_controls().energy_tolerance_j);
            close(
                report.absolute_energy_residual_j[j],
                0.0,
                piston_controls().energy_tolerance_j,
            );
        }
        close(middle.heat.energy_residual_j, 0.0, 3e-12);
        close(
            report.total_energy_residual_j,
            0.0,
            piston_controls().energy_tolerance_j,
        );
        close(
            report.combined_momentum_residual_n_s,
            0.0,
            piston_controls().momentum_tolerance_n_s,
        );
        drive += report.piston.external_force_work_j;
    }
    close(
        physical_change(left.flow(), before[0])
            + physical_change(right.flow(), before[1])
            + WALL_CAPACITY * (wall_temperature(&thermal) - 4.0)
            + 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2)
            - initial_kinetic
            - drive,
        0.0,
        3e-7,
    );
    assert!(left.flow().mesh().volumes().iter().sum::<f64>() > 0.8 * AREA);
    assert!(right.flow().mesh().volumes().iter().sum::<f64>() < 1.2 * AREA);
    assert!(left.flow().totals().species_masses[1] > 0.0);
    assert!(right.flow().totals().species_masses[1] > 0.0);
    assert!(wall_temperature(&thermal) < 4.0);
}

#[test]
fn g1_shortened_thermal_interval_replays_both_chemistries_at_moved_volumes() {
    let make = || {
        (
            chamber(&mechanism(40.0, 500.0, 0.0), 0.0, 1.0, 4, chemistry()),
            chamber(&mechanism(0.0, 2.0, 30.0), 1.0, 1.0, 4, chemistry()),
        )
    };
    let (mut left, mut right) = make();
    let mut thermal = wall(&left, &right, 10.0, [3.0, 1.0]);
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let report = advance(
        &mut left,
        &mut right,
        &mut thermal,
        &mut piston,
        0.2,
        || false,
    )
    .unwrap();
    assert!(report.coupling_retries > 0);
    let dt = report.elapsed.value();

    let (mut replay_left, mut replay_right) = make();
    let mut replay_wall = wall(&replay_left, &replay_right, 10.0, [3.0, 1.0]);
    let volumes = [
        replay_left.flow().mesh().volumes().to_vec(),
        replay_right.flow().mesh().volumes().to_vec(),
    ];
    let mut source = [ChemistryReport::default(); 2];
    for (j, flow) in [&mut replay_left, &mut replay_right]
        .into_iter()
        .enumerate()
    {
        react_cells(
            &mut flow.accepted,
            &mut flow.reactor,
            &mut flow.amounts,
            0.5 * dt,
            &mut source[j],
            &mut || false,
        )
        .unwrap();
    }
    let mut replay_piston = Piston::new(0.2, 0.0).unwrap();
    let middle: PistonStep<SharedPistonThermalStepReport> = replay_left
        .accepted
        .step_with_shared_piston_and_heat(
            &mut replay_right.accepted,
            &mut replay_wall,
            &mut replay_piston,
            Time::new(dt),
            StepControl::default(),
            heat_controls(),
            0.01,
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.transport.left.elapsed, report.elapsed);
    for (j, flow) in [&mut replay_left, &mut replay_right]
        .into_iter()
        .enumerate()
    {
        assert_ne!(flow.flow().mesh().volumes(), volumes[j]);
        react_cells(
            &mut flow.accepted,
            &mut flow.reactor,
            &mut flow.amounts,
            0.5 * dt,
            &mut source[j],
            &mut || false,
        )
        .unwrap();
        for k in 0..2 {
            close(
                report.chemistry[j].species_change_kg[k],
                source[j].species_change_kg[k],
                2e-10,
            );
        }
    }
    compare_fields(&left, &replay_left, 3e-8);
    compare_fields(&right, &replay_right, 3e-8);
    close(piston.velocity_m_s(), replay_piston.velocity_m_s(), 2e-9);
    close(
        wall_temperature(&thermal),
        wall_temperature(&replay_wall),
        2e-9,
    );
}

#[test]
fn g1_chemical_gauges_and_heat_preserve_independent_isomerization_amounts() {
    let make = |references: [f64; 2]| {
        (
            chamber(
                &mechanism(0.0, 5.0, references[0]),
                0.0,
                0.8,
                3,
                chemistry(),
            ),
            chamber(
                &mechanism(0.0, 2.0, references[1]),
                0.8,
                1.2,
                4,
                chemistry(),
            ),
        )
    };
    let (mut left, mut right) = make([0.0, 0.0]);
    let (mut shifted_left, mut shifted_right) = make([20_000.0, -30_000.0]);
    let mut thermal = wall(&left, &right, 4.0, [0.9, 0.3]);
    let mut shifted_wall = wall(&shifted_left, &shifted_right, 4.0, [0.9, 0.3]);
    let mut piston = Piston::new(0.2, 0.04).unwrap();
    let mut shifted_piston = piston;
    let masses = [left.flow().totals().mass(), right.flow().totals().mass()];
    for _ in 0..8 {
        advance(
            &mut left,
            &mut right,
            &mut thermal,
            &mut piston,
            0.002,
            || false,
        )
        .unwrap();
        advance(
            &mut shifted_left,
            &mut shifted_right,
            &mut shifted_wall,
            &mut shifted_piston,
            0.002,
            || false,
        )
        .unwrap();
    }
    // Equal original calorics make A <-> B's exact amount trajectory independent
    // of temperature, volume and pressure: n_B/n = (1-exp(-2*k*t))/2.
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        let time = flow.flow().time().value();
        let expected = 0.5 * masses[j] * (1.0 - (-2.0 * [5.0, 2.0][j] * time).exp());
        close(flow.flow().totals().species_masses[1], expected, 3e-7);
    }
    compare_fields(&left, &shifted_left, 4e-8);
    compare_fields(&right, &shifted_right, 4e-8);
    close(piston.velocity_m_s(), shifted_piston.velocity_m_s(), 3e-9);
    close(
        wall_temperature(&thermal),
        wall_temperature(&shifted_wall),
        3e-9,
    );
}

#[test]
fn g1_suspended_reacting_piston_returns_only_the_authored_damper_fraction() {
    for fraction in [0.0, 0.35, 1.0] {
        let (mut left, mut right) = pair();
        let mut thermal = wall(&left, &right, 1.4, [3.0, 0.2]);
        let mut piston = Piston::new(0.4, 0.2).unwrap();
        let suspension = PistonSuspension::new(8.0, 2.5, 0.75).unwrap();
        let before = [left.flow().totals(), right.flow().totals()];
        let initial_kinetic = 0.5 * 0.4 * 0.2_f64.powi(2);
        let initial_spring = suspension.stored_energy_j(0.8).unwrap();
        let (mut loss, mut drive) = (0.0, 0.0);
        for _ in 0..8 {
            let x0 = left.flow().mesh().positions()[left.flow().cell_count()];
            let v0 = piston.velocity_m_s();
            let report = left
                .step_with_shared_suspended_piston_and_heat(
                    &mut right,
                    &mut thermal,
                    &mut piston,
                    Time::new(0.002),
                    ReactiveStepControl::default(),
                    heat_controls(),
                    0.02,
                    suspension,
                    fraction,
                    piston_controls(),
                    || false,
                )
                .unwrap();
            let mechanical = report.piston;
            let dt = report.elapsed.value();
            let midpoint = 0.5 * (v0 + piston.velocity_m_s());
            let expected_midpoint = (2.0 * 0.4 * v0
                + mechanical.piston.gas_impulse_n_s
                + dt * (0.02 - 8.0 * (x0 - 0.75)))
                / (2.0 * 0.4 + 2.5 * dt + 0.5 * 8.0 * dt * dt);
            close(midpoint, expected_midpoint, 4e-11);
            let damping = 2.5 * midpoint * midpoint * dt;
            let heat = mechanical.piston.gas.heat;
            close(mechanical.damping_energy_j, damping, 3e-14);
            close(heat.damping_heat_to_wall_j, fraction * damping, 3e-14);
            close(
                heat.damping_loss_external_j,
                (1.0 - fraction) * damping,
                3e-14,
            );
            close(heat.energy_residual_j, 0.0, 3e-12);
            close(
                report.total_energy_residual_j,
                0.0,
                piston_controls().energy_tolerance_j,
            );
            loss += heat.damping_loss_external_j;
            drive += mechanical.applied_force_work_j;
        }
        let endpoint = left.flow().mesh().positions()[left.flow().cell_count()];
        close(
            physical_change(left.flow(), before[0])
                + physical_change(right.flow(), before[1])
                + WALL_CAPACITY * (wall_temperature(&thermal) - 1.4)
                + 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2)
                - initial_kinetic
                + suspension.stored_energy_j(endpoint).unwrap()
                - initial_spring
                + loss
                - drive,
            0.0,
            3e-7,
        );
    }
}

#[test]
fn g4_final_chemistry_refusal_cannot_publish_wall_heat_or_damper_deposit() {
    let source = mechanism(0.0, 2.0, 0.0);
    let narrow = IsochoricKineticsSolve {
        maximum_temperature_k: 1.000001,
        ..chemistry()
    };
    let mut left = chamber(&source, 0.0, 1.0, 8, chemistry());
    let mut right = chamber(&source, 1.0, 1.0, 8, narrow);
    let mut thermal = wall(&left, &right, 4.0, [0.5, 0.5]);
    let mut piston = Piston::new(0.4, 0.04).unwrap();
    let original_piston = piston;
    let original_wall = thermal.storage().thermal;
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let error = left
        .step_with_shared_suspended_piston_and_heat(
            &mut right,
            &mut thermal,
            &mut piston,
            Time::new(0.002),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            heat_controls(),
            0.01,
            PistonSuspension::new(2.0, 0.8, 0.9).unwrap(),
            1.0,
            piston_controls(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(
        error, ReactiveError::RetryLimit { last, .. }
        if *last == ReactiveError::Invalid("cell temperature is outside chemistry interval")
    ));
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.accepted.states, states[j]);
        assert_eq!(flow.flow().mesh(), &meshes[j]);
        assert_eq!(flow.flow().time(), Time::new(0.0));
    }
    assert_eq!(piston, original_piston);
    assert_eq!(thermal.storage().thermal, original_wall);
    assert_eq!(thermal.last_report(), None);
}

#[test]
fn g4_cancellation_after_prior_thermal_motion_keeps_all_owners_and_retry() {
    let prepare = || {
        let (mut left, mut right) = pair();
        let mut thermal = wall(&left, &right, 4.0, [0.9, 0.3]);
        let mut piston = Piston::new(0.2, 0.04).unwrap();
        advance(
            &mut left,
            &mut right,
            &mut thermal,
            &mut piston,
            0.002,
            || false,
        )
        .unwrap();
        (left, right, thermal, piston)
    };
    let (mut expected_left, mut expected_right, mut expected_wall, mut expected_piston) = prepare();
    let states = [
        expected_left.accepted.states.clone(),
        expected_right.accepted.states.clone(),
    ];
    let meshes = [
        expected_left.flow().mesh().clone(),
        expected_right.flow().mesh().clone(),
    ];
    let time = expected_left.flow().time();
    let original_piston = expected_piston;
    let original_wall = expected_wall.storage().thermal;
    let last_heat = expected_wall.last_report();
    let mut polls = 0;
    let expected = advance(
        &mut expected_left,
        &mut expected_right,
        &mut expected_wall,
        &mut expected_piston,
        0.002,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let (mut left, mut right, mut thermal, mut piston) = prepare();
        let mut seen = 0;
        assert_eq!(
            advance(
                &mut left,
                &mut right,
                &mut thermal,
                &mut piston,
                0.002,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(ReactiveError::Cancelled),
        );
        assert_eq!(seen, cutoff);
        for (j, flow) in [&left, &right].into_iter().enumerate() {
            assert_eq!(flow.accepted.states, states[j]);
            assert_eq!(flow.flow().mesh(), &meshes[j]);
            assert_eq!(flow.flow().time(), time);
        }
        assert_eq!(piston, original_piston);
        assert_eq!(thermal.storage().thermal, original_wall);
        assert_eq!(thermal.last_report(), last_heat);
        let recovered = advance(
            &mut left,
            &mut right,
            &mut thermal,
            &mut piston,
            0.002,
            || false,
        )
        .unwrap();
        assert_eq!(recovered.elapsed, expected.elapsed);
        compare_fields(&left, &expected_left, 3e-8);
        compare_fields(&right, &expected_right, 3e-8);
        close(piston.velocity_m_s(), expected_piston.velocity_m_s(), 2e-9);
        close(
            wall_temperature(&thermal),
            wall_temperature(&expected_wall),
            2e-9,
        );
    }
}
