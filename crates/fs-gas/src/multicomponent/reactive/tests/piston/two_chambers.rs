use super::*;
use crate::multicomponent::Mesh1d;

mod molecular;

fn chamber(
    mechanism: &IdealGasKineticsV1,
    start: f64,
    length: f64,
    cells: usize,
    chemistry: IsochoricKineticsSolve,
) -> ReactingEuler1d {
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let positions: Vec<_> = (0..=cells)
        .map(|i| start + length * (i as f64 / cells as f64))
        .collect();
    let volumes = positions.windows(2).map(|p| AREA * (p[1] - p[0])).collect();
    let mesh = Mesh1d::new(positions, vec![AREA; cells + 1], volumes).unwrap();
    let state = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let flow = MixtureEuler1d::new_on_mesh(
        gas,
        mesh,
        vec![state; cells],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    ReactingEuler1d::new(flow, mechanism.clone(), chemistry, || false).unwrap()
}

fn pair() -> (ReactingEuler1d, ReactingEuler1d) {
    (
        chamber(&mechanism(8.0, 30.0, 100.0), 0.0, 0.8, 8, chemistry()),
        chamber(&mechanism(2.0, 8.0, -30.0), 0.8, 1.2, 12, chemistry()),
    )
}

fn step(
    left: &mut ReactingEuler1d,
    right: &mut ReactingEuler1d,
    piston: &mut Piston,
    duration: Time,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveSharedPistonStep, ReactiveError> {
    left.step_with_shared_piston(
        right,
        piston,
        duration,
        ReactiveStepControl::default(),
        0.0,
        piston_controls(),
        cancelled,
    )
}

fn physical_change(flow: &MixtureEuler1d, before: Totals) -> f64 {
    let after = flow.totals();
    sum(
        core::iter::once(after.energy - before.energy).chain((0..flow.gas().species_count()).map(
            |k| {
                (after.species_masses[k] - before.species_masses[k])
                    * flow.gas().energy_references()[k]
            },
        )),
    )
}

#[test]
fn g1_independent_reacting_chambers_drive_one_piston_and_close_physical_energy() {
    let (mut left, mut right) = pair();
    let before = [left.flow().totals(), right.flow().totals()];
    let mechanisms = [left.mechanism().clone(), right.mechanism().clone()];
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let mut support = 0.0;
    let mut work = 0.0;
    let mut imposed_impulse = 0.0;
    let mut sources = [[0.0; 2]; 2];
    let force = -0.001;
    for _ in 0..30 {
        let report = left
            .step_with_shared_piston(
                &mut right,
                &mut piston,
                Time::new(0.002),
                ReactiveStepControl::default(),
                force,
                piston_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report.left.elapsed, report.right.elapsed);
        assert_eq!(report.left.elapsed, report.piston.gas.left.elapsed);
        assert_eq!(report.right.elapsed, report.piston.gas.right.elapsed);
        assert_eq!(left.flow().time(), right.flow().time());
        for (j, side) in [&report.left, &report.right].into_iter().enumerate() {
            assert_eq!(side.chemistry.cell_solves, [16, 24][j]);
            assert!(side.chemistry.referenced_energy_change_j > 0.0);
            close(side.mass_residual_kg, 0.0, 3e-12);
            for k in 0..2 {
                assert_eq!(side.transport.boundary_change.species_masses[k], 0.0);
                close(side.conservation_residual.species_masses[k], 0.0, 3e-12);
                sources[j][k] += side.chemistry.species_change_kg[k];
            }
        }
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
        assert_eq!(report.piston.exterior_pressure_work_j, 0.0);
        support += report.piston.fixed_wall_impulse_n_s;
        work += report.piston.external_force_work_j;
        imposed_impulse += force * report.left.elapsed.value();
    }
    let after = [left.flow().totals(), right.flow().totals()];
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    close(
        physical_change(left.flow(), before[0])
            + physical_change(right.flow(), before[1])
            + kinetic
            - work,
        0.0,
        3e-7,
    );
    close(
        after[0].momentum - before[0].momentum + after[1].momentum - before[1].momentum
            + piston.mass_kg() * piston.velocity_m_s()
            + support
            - imposed_impulse,
        0.0,
        3e-10,
    );
    for (j, flow) in [&left, &right].into_iter().enumerate() {
        assert_eq!(flow.mechanism(), &mechanisms[j]);
        for k in 0..2 {
            close(
                after[j].species_masses[k] - before[j].species_masses[k],
                sources[j][k],
                3e-12,
            );
        }
        assert!(after[j].species_masses[1] > 0.0);
        assert!(after[j].energy > before[j].energy);
    }
    assert!(
        piston.velocity_m_s() > 0.001,
        "stronger left reaction must drive the shared piston"
    );
    let interface = left.flow().mesh().positions()[left.flow().cell_count()];
    assert!(interface > 0.8);
    assert_eq!(interface, right.flow().mesh().positions()[0]);
    assert_eq!(left.flow().mesh().positions()[0], 0.0);
    assert_eq!(
        right.flow().mesh().positions()[right.flow().cell_count()],
        2.0
    );
    close(
        left.flow().mesh().volumes().iter().sum::<f64>()
            + right.flow().mesh().volumes().iter().sum::<f64>(),
        2.0 * AREA,
        3e-14,
    );
}

#[test]
fn g1_heat_release_restarts_both_sources_at_common_duration_and_actual_moved_volumes() {
    let make = || {
        (
            chamber(&mechanism(40.0, 500.0, 0.0), 0.0, 1.0, 4, chemistry()),
            chamber(&mechanism(0.0, 2.0, 30.0), 1.0, 1.0, 4, chemistry()),
        )
    };
    let (mut left, mut right) = make();
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let report = step(&mut left, &mut right, &mut piston, Time::new(0.2), || false).unwrap();
    assert!(
        report.left.coupling_retries > 0,
        "chemistry must shorten the original acoustic bound"
    );
    assert_eq!(report.left.coupling_retries, report.right.coupling_retries);
    let dt = report.left.elapsed.value();
    let (mut replay_left, mut replay_right) = make();
    let volumes = [
        replay_left.flow().mesh().volumes().to_vec(),
        replay_right.flow().mesh().volumes().to_vec(),
    ];
    let mut chemistry = [ChemistryReport::default(); 2];
    for (j, flow) in [&mut replay_left, &mut replay_right]
        .into_iter()
        .enumerate()
    {
        react_cells(
            &mut flow.accepted,
            &mut flow.reactor,
            &mut flow.amounts,
            0.5 * dt,
            &mut chemistry[j],
            &mut || false,
        )
        .unwrap();
    }
    let mut replay_piston = Piston::new(0.2, 0.0).unwrap();
    let middle = replay_left
        .accepted
        .step_with_shared_piston(
            &mut replay_right.accepted,
            &mut replay_piston,
            Time::new(dt),
            StepControl::default(),
            0.0,
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.left.elapsed, report.left.elapsed);
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
            &mut chemistry[j],
            &mut || false,
        )
        .unwrap();
    }
    close(piston.velocity_m_s(), replay_piston.velocity_m_s(), 2e-9);
    for (j, (actual, expected)) in [(&left, &replay_left), (&right, &replay_right)]
        .into_iter()
        .enumerate()
    {
        assert_eq!(actual.flow().time(), expected.flow().time());
        for (&a, &b) in actual.accepted.states.iter().zip(&expected.accepted.states) {
            close(a, b, 2e-9);
        }
        for (&a, &b) in actual
            .flow()
            .mesh()
            .volumes()
            .iter()
            .zip(expected.flow().mesh().volumes())
        {
            close(a, b, 2e-12);
        }
        let source = [&report.left, &report.right][j];
        for k in 0..2 {
            close(
                source.chemistry.species_change_kg[k],
                chemistry[j].species_change_kg[k],
                2e-10,
            );
        }
    }
}

#[test]
fn g4_right_postcompression_source_refusal_preserves_both_accepted_chambers() {
    let source = mechanism(0.0, 2.0, 0.0);
    let narrow = IsochoricKineticsSolve {
        maximum_temperature_k: 1.0000001,
        ..chemistry()
    };
    let make = || {
        (
            chamber(&source, 0.0, 1.0, 8, chemistry()),
            chamber(&source, 1.0, 1.0, 8, narrow),
        )
    };
    let dt = Time::new(0.001);
    let (mut stages_left, mut stages_right) = make();
    for flow in [&mut stages_left, &mut stages_right] {
        react_cells(
            &mut flow.accepted,
            &mut flow.reactor,
            &mut flow.amounts,
            0.5 * dt.value(),
            &mut ChemistryReport::default(),
            &mut || false,
        )
        .unwrap();
    }
    let mut trial_piston = Piston::new(0.1, 0.0).unwrap();
    let middle = stages_left
        .accepted
        .step_with_shared_piston(
            &mut stages_right.accepted,
            &mut trial_piston,
            dt,
            StepControl::default(),
            1.0,
            piston_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(middle.gas.left.elapsed, dt);
    assert!(middle.displacement_m > 0.0);
    react_cells(
        &mut stages_left.accepted,
        &mut stages_left.reactor,
        &mut stages_left.amounts,
        0.5 * dt.value(),
        &mut ChemistryReport::default(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        react_cells(
            &mut stages_right.accepted,
            &mut stages_right.reactor,
            &mut stages_right.amounts,
            0.5 * dt.value(),
            &mut ChemistryReport::default(),
            &mut || false,
        ),
        Err(ReactiveError::Invalid(
            "cell temperature is outside chemistry interval"
        )),
    );
    let (mut left, mut right) = make();
    let states = [left.accepted.states.clone(), right.accepted.states.clone()];
    let meshes = [left.flow().mesh().clone(), right.flow().mesh().clone()];
    let mut piston = Piston::new(0.1, 0.0).unwrap();
    let original = piston;
    let error = left
        .step_with_shared_piston(
            &mut right,
            &mut piston,
            dt,
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
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
    assert_eq!(piston, original);
}

#[test]
fn g4_final_cancel_after_prior_motion_retains_both_sources_meshes_and_retry() {
    let dt = Time::new(0.002);
    let (mut reference_left, mut reference_right) = pair();
    let mut reference_piston = Piston::new(0.2, 0.0).unwrap();
    step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_piston,
        dt,
        || false,
    )
    .unwrap();
    let states = [
        reference_left.accepted.states.clone(),
        reference_right.accepted.states.clone(),
    ];
    let meshes = [
        reference_left.flow().mesh().clone(),
        reference_right.flow().mesh().clone(),
    ];
    let time = reference_left.flow().time();
    let original_piston = reference_piston;
    let mut polls = 0;
    let expected = step(
        &mut reference_left,
        &mut reference_right,
        &mut reference_piston,
        dt,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let (mut left, mut right) = pair();
        let mut piston = Piston::new(0.2, 0.0).unwrap();
        step(&mut left, &mut right, &mut piston, dt, || false).unwrap();
        let mut seen = 0;
        assert_eq!(
            step(&mut left, &mut right, &mut piston, dt, || {
                seen += 1;
                seen == cutoff
            }),
            Err(ReactiveError::Cancelled),
        );
        assert_eq!(seen, cutoff);
        for (j, flow) in [&left, &right].into_iter().enumerate() {
            assert_eq!(flow.accepted.states, states[j]);
            assert_eq!(flow.flow().mesh(), &meshes[j]);
            assert_eq!(flow.flow().time(), time);
        }
        assert_eq!(piston, original_piston);
        let recovered = step(&mut left, &mut right, &mut piston, dt, || false).unwrap();
        assert_eq!(recovered.left.elapsed, expected.left.elapsed);
        close(
            piston.velocity_m_s(),
            reference_piston.velocity_m_s(),
            2e-10,
        );
        for (flow, reference) in [(&left, &reference_left), (&right, &reference_right)] {
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
}
