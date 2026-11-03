use super::*;
use crate::dynamics::piston::{Piston, PistonControl, PistonLoad};
use crate::multicomponent::Boundary;

mod molecular;
mod thermal;
mod inertial_walls;
mod two_chambers;

const AREA: f64 = 0.1;

fn closed(
    mechanism: &IdealGasKineticsV1,
    cells: usize,
    chemistry: IsochoricKineticsSolve,
) -> ReactingEuler1d {
    let mut flow = field(mechanism, cells, 0.0, |_| [1.0, 0.0]);
    flow.boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    ReactingEuler1d::new(flow, mechanism.clone(), chemistry, || false).unwrap()
}

fn piston_controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-1.0, 1.0],
        momentum_tolerance_n_s: 2e-11,
        energy_tolerance_j: 2e-8,
        maximum_iterations: 64,
    }
}

fn load() -> PistonLoad {
    PistonLoad {
        exterior_pressure_pa: 1.0,
        external_force_n: 0.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_chemical_energy_drives_resolved_pressure_and_piston_without_added_heat() {
    let mechanism = mechanism(8.0, 30.0, 0.0);
    let mut flow = closed(&mechanism, 16, chemistry());
    let before = flow.flow().totals();
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let mut exterior_work = 0.0;
    let mut exterior_impulse = 0.0;
    let mut support = 0.0;
    let mut species_change = [0.0; 2];
    let mut steps = 0;
    while flow.flow().time().value() < 0.06 {
        let report = flow
            .step_with_piston(
                &mut piston,
                Time::new(0.06 - flow.flow().time().value()),
                ReactiveStepControl::default(),
                load(),
                piston_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report.flow.elapsed, report.piston.gas.elapsed);
        assert_eq!(report.flow.chemistry.cell_solves, 32);
        assert!(report.flow.chemistry.referenced_energy_change_j > 0.0);
        assert!(report.piston.displacement_m > 0.0);
        assert!(report.piston.gas.right_boundary_change.energy < 0.0);
        assert_eq!(report.piston.gas.left_boundary_change.energy, 0.0);
        close(report.flow.mass_residual_kg, 0.0, 3e-12);
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
        for (k, change) in species_change.iter_mut().enumerate() {
            assert_eq!(report.flow.transport.boundary_change.species_masses[k], 0.0);
            close(
                report.flow.conservation_residual.species_masses[k],
                0.0,
                3e-12,
            );
            *change += report.flow.chemistry.species_change_kg[k];
        }
        exterior_work += report.piston.exterior_pressure_work_j;
        exterior_impulse -= AREA * report.flow.elapsed.value();
        support += report.piston.fixed_wall_impulse_n_s;
        steps += 1;
        assert!(steps < 500);
    }
    let after = flow.flow().totals();
    let reference_change: f64 = (0..2)
        .map(|k| {
            (after.species_masses[k] - before.species_masses[k])
                * flow.flow().gas().energy_references()[k]
        })
        .sum();
    let physical_change = (after.energy - before.energy) + reference_change;
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    close(physical_change + kinetic + exterior_work, 0.0, 3e-8);
    close(
        after.momentum - before.momentum + piston.mass_kg() * piston.velocity_m_s(),
        exterior_impulse - support,
        3e-10,
    );
    close(after.mass(), before.mass(), 3e-12);
    for (k, &change) in species_change.iter().enumerate() {
        close(
            after.species_masses[k] - before.species_masses[k],
            change,
            3e-12,
        );
    }
    assert!(
        reference_change < -0.1,
        "chemical storage must supply the energy"
    );
    assert!(
        after.energy > before.energy + 0.1,
        "sensible energy must rise without a second heat source"
    );
    assert!(after.species_masses[0] < before.species_masses[0]);
    assert!(after.species_masses[1] > 0.0);
    assert!(piston.velocity_m_s() > 0.01);
    assert_eq!(flow.flow().mesh().positions()[0], 0.0);
    assert!(flow.flow().mesh().positions()[16] > 1.0);
    close(
        flow.flow().mesh().volumes().iter().sum(),
        AREA * flow.flow().mesh().positions()[16],
        5e-14,
    );
    assert!(
        flow.flow()
            .gas()
            .primitive(flow.flow().cell(0).unwrap(), || false)
            .unwrap()
            .pressure
            > 1.1
    );
}

#[test]
fn g1_acoustic_restart_replays_both_sources_at_accepted_dt_and_moved_volume() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = closed(&mechanism, 4, chemistry());
    let mut piston = Piston::new(2.0, 0.0).unwrap();
    let report = flow
        .step_with_piston(
            &mut piston,
            Time::new(0.2),
            ReactiveStepControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
    assert!(
        report.flow.coupling_retries > 0,
        "heat release must shorten the common acoustic duration"
    );
    assert_eq!(report.flow.elapsed, report.piston.gas.elapsed);
    assert_eq!(flow.flow().time(), report.flow.elapsed);
    assert_eq!(report.flow.chemistry.cell_solves, 8);
    let dt = report.flow.elapsed.value();
    // Replay with the original public physical stages at the final duration.
    // In particular, C2 receives the actual moved volumes, not the original V.
    let mut replay = closed(&mechanism, 4, chemistry());
    let original_volume = replay.flow().mesh().volumes().to_vec();
    let mut wall = Piston::new(2.0, 0.0).unwrap();
    let mut sources = ChemistryReport::default();
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut sources,
        &mut || false,
    )
    .unwrap();
    let middle = replay
        .accepted
        .step_with_piston(
            &mut wall,
            Time::new(dt),
            StepControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.elapsed, report.flow.elapsed);
    assert_ne!(replay.flow().mesh().volumes(), original_volume);
    react_cells(
        &mut replay.accepted,
        &mut replay.reactor,
        &mut replay.amounts,
        0.5 * dt,
        &mut sources,
        &mut || false,
    )
    .unwrap();
    close(piston.velocity_m_s(), wall.velocity_m_s(), 2e-9);
    for (&a, &b) in flow.accepted.states.iter().zip(&replay.accepted.states) {
        close(a, b, 2e-9);
    }
    for (&a, &b) in flow
        .flow()
        .mesh()
        .volumes()
        .iter()
        .zip(replay.flow().mesh().volumes())
    {
        close(a, b, 2e-12);
    }
    for k in 0..2 {
        close(
            report.flow.chemistry.species_change_kg[k],
            sources.species_change_kg[k],
            2e-10,
        );
    }
    close(
        report.total_energy_residual_j,
        0.0,
        piston_controls().energy_tolerance_j,
    );
}

#[test]
fn g4_second_chemistry_refusal_rolls_back_successful_pressure_and_mesh_motion() {
    let mechanism = mechanism(0.0, 2.0, 0.0);
    let chemistry = IsochoricKineticsSolve {
        maximum_temperature_k: 1.0000001,
        ..chemistry()
    };
    let duration = Time::new(0.001);
    let load = PistonLoad {
        exterior_pressure_pa: 2.0,
        external_force_n: 0.0,
    };
    // First prove that the original first source and pressure step are admitted,
    // while adiabatic compression leaves the narrower reactor interval.
    let mut stages = closed(&mechanism, 8, chemistry);
    let mut stage_piston = Piston::new(0.1, 0.0).unwrap();
    let mut sources = ChemistryReport::default();
    react_cells(
        &mut stages.accepted,
        &mut stages.reactor,
        &mut stages.amounts,
        0.5 * duration.value(),
        &mut sources,
        &mut || false,
    )
    .unwrap();
    let middle = stages
        .accepted
        .step_with_piston(
            &mut stage_piston,
            duration,
            StepControl::default(),
            load,
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.elapsed, duration);
    assert!(middle.displacement_m < 0.0);
    assert_eq!(
        react_cells(
            &mut stages.accepted,
            &mut stages.reactor,
            &mut stages.amounts,
            0.5 * duration.value(),
            &mut sources,
            &mut || false
        ),
        Err(ReactiveError::Invalid(
            "cell temperature is outside chemistry interval"
        ))
    );
    let mut flow = closed(&mechanism, 8, chemistry);
    let before = flow.accepted.states.clone();
    let mesh = flow.flow().mesh().clone();
    let mut piston = Piston::new(0.1, 0.0).unwrap();
    let original_piston = piston;
    let result = flow.step_with_piston(
        &mut piston,
        duration,
        ReactiveStepControl {
            maximum_coupling_retries: 0,
            ..ReactiveStepControl::default()
        },
        load,
        piston_controls(),
        || false,
    );
    assert!(matches!(result, Err(ReactiveError::RetryLimit { last, .. })
        if *last == ReactiveError::Invalid("cell temperature is outside chemistry interval")));
    assert_eq!(flow.accepted.states, before);
    assert_eq!(flow.flow().mesh(), &mesh);
    assert_eq!(flow.flow().time(), Time::new(0.0));
    assert_eq!(piston, original_piston);
}

#[test]
fn g4_cancellation_after_prior_motion_preserves_every_owner_and_reuses_scratch() {
    let mechanism = mechanism(10.0, 30.0, 0.0);
    let duration = Time::new(0.002);
    let mut reference = closed(&mechanism, 8, chemistry());
    let mut reference_piston = Piston::new(0.2, 0.0).unwrap();
    reference
        .step_with_piston(
            &mut reference_piston,
            duration,
            ReactiveStepControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
    let before = reference.accepted.states.clone();
    let mesh = reference.flow().mesh().clone();
    let time = reference.flow().time();
    let original_piston = reference_piston;
    let mut polls = 0;
    let expected = reference
        .step_with_piston(
            &mut reference_piston,
            duration,
            ReactiveStepControl::default(),
            load(),
            piston_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, 8, polls / 2, polls] {
        let mut flow = closed(&mechanism, 8, chemistry());
        let mut piston = Piston::new(0.2, 0.0).unwrap();
        flow.step_with_piston(
            &mut piston,
            duration,
            ReactiveStepControl::default(),
            load(),
            piston_controls(),
            || false,
        )
        .unwrap();
        let mut seen = 0;
        assert_eq!(
            flow.step_with_piston(
                &mut piston,
                duration,
                ReactiveStepControl::default(),
                load(),
                piston_controls(),
                || {
                    seen += 1;
                    seen == cutoff
                }
            ),
            Err(ReactiveError::Cancelled)
        );
        assert_eq!(seen, cutoff);
        assert_eq!(flow.accepted.states, before);
        assert_eq!(flow.flow().mesh(), &mesh);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(piston, original_piston);
        let recovered = flow
            .step_with_piston(
                &mut piston,
                duration,
                ReactiveStepControl::default(),
                load(),
                piston_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(recovered.flow.elapsed, expected.flow.elapsed);
        close(
            piston.velocity_m_s(),
            reference_piston.velocity_m_s(),
            2e-10,
        );
        for (&a, &b) in flow.accepted.states.iter().zip(&reference.accepted.states) {
            close(a, b, 2e-10);
        }
        for (&a, &b) in flow
            .flow()
            .mesh()
            .positions()
            .iter()
            .zip(reference.flow().mesh().positions())
        {
            close(a, b, 2e-13);
        }
    }
}
